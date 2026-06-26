"""Gaussian likelihood over an arbitrary collection of datasets.

The user composes a :class:`GaussianLikelihood` from whatever
:class:`~halogal.inference.data.Dataset` objects they want to fit -- any mix of
UVLF, clustering, and number-density terms, across one or several redshift bins.
That composition *is* the choice of what data the fit constrains.

For efficiency the likelihood keeps one :class:`~halogal.model.Observables` per
redshift bin so that the expensive ``halomod`` correlation model is initialised
once and only *updated* during sampling. These live objects are deliberately
excluded from pickling (:meth:`__getstate__`): under multiprocessing/MPI each
worker reconstructs its own, because ``halomod`` state is not safely picklable.
"""

from collections import OrderedDict, defaultdict

import numpy as np

__all__ = ["GaussianLikelihood"]


class GaussianLikelihood:
    """Total Gaussian log-likelihood for a set of datasets.

    Parameters
    ----------
    datasets : list of Dataset
        The observables to fit. At most one clustering dataset is allowed per
        redshift bin (a single cached ``halomod`` model serves each bin).
    parametrization : ModelParametrization
        Maps a named parameter dict to per-redshift model parameters.
    catch_errors : bool, optional
        If True (default), exceptions raised while evaluating a model (e.g. an
        occasional ``halomod`` failure at extreme parameters) are caught and the
        likelihood returns ``-inf`` for that point. Set False to surface errors
        during development.
    """

    def __init__(self, datasets, parametrization, catch_errors=True):
        self.datasets = list(datasets)
        if not self.datasets:
            raise ValueError("GaussianLikelihood requires at least one Dataset")
        self.parametrization = parametrization
        self.catch_errors = catch_errors
        self._validate()
        # Lazily constructed per process; see _ensure_setup / __getstate__.
        self._bins = None

    def _validate(self):
        clustering_per_z = defaultdict(int)
        for d in self.datasets:
            if d.kind == "clustering":
                clustering_per_z[d.z] += 1
        for z, count in clustering_per_z.items():
            if count > 1:
                raise ValueError(
                    f"At most one clustering dataset per redshift bin "
                    f"(z={z:g} has {count}); a single halomod model is cached "
                    f"per bin."
                )

    @property
    def redshifts(self):
        """Sorted unique redshifts spanned by the datasets."""
        return sorted({d.z for d in self.datasets})

    def _ensure_setup(self):
        """Build the per-bin Observables (idempotent, per process)."""
        if self._bins is not None:
            return
        from ..model import HODModel, Observables

        bins = OrderedDict()
        for d in self.datasets:
            if d.z not in bins:
                obs = Observables(HODModel(z=d.z))
                bins[d.z] = {"obs": obs, "datasets": [], "clustering": None}
            bins[d.z]["datasets"].append(d)
            if d.kind == "clustering":
                bins[d.z]["clustering"] = d

        # Initialise each bin's correlation model once (expensive).
        for z, b in bins.items():
            d = b["clustering"]
            if d is not None:
                b["obs"].initialize_correlation_model(
                    MUV_thresh1=d.MUV_thresh1,
                    MUV_thresh2=d.MUV_thresh2,
                    **d.corr_kwargs,
                )
        self._bins = bins

    def _accumulate(self, named_params, reduce_fn):
        """Shared loop over bins/datasets applying ``reduce_fn(dataset, model)``.

        Returns the summed reductions, or ``-inf`` if any model is non-finite
        (when reducing log-likelihoods) / on a caught error.
        """
        self._ensure_setup()
        total = 0.0
        for z, b in self._bins.items():
            model_params = self.parametrization.params_for_redshift(
                named_params, z
            )
            for d in b["datasets"]:
                model = d.predict(b["obs"], model_params)
                total += reduce_fn(d, model)
        return total

    def log_likelihood(self, named_params):
        """Total Gaussian log-likelihood for a named parameter dict.

        Parameters
        ----------
        named_params : dict
            Sampling-space parameter values (e.g. from
            ``ParameterSet.to_dict(theta)``).

        Returns
        -------
        float
            Sum of per-dataset Gaussian log-likelihoods, or ``-inf`` if any
            prediction is non-finite or (optionally) raises.
        """
        def reduce_fn(d, model):
            ll = d.log_likelihood(model)
            if not np.isfinite(ll):
                raise _NonFinite()
            return ll

        try:
            return self._accumulate(named_params, reduce_fn)
        except _NonFinite:
            return -np.inf
        except Exception:
            if self.catch_errors:
                return -np.inf
            raise

    def chi2(self, named_params):
        """Total generalised chi-square (diagnostic, no normalisation)."""
        return self._accumulate(
            named_params, lambda d, model: d.chi2(model)
        )

    def predict_all(self, named_params):
        """Model predictions for every dataset, keyed by label.

        Useful for overlaying a best-fit model on the data.

        Returns
        -------
        dict
            ``{label: {"x": abscissa, "data": measured, "model": predicted,
            "z": redshift, "kind": kind}}``.
        """
        self._ensure_setup()
        out = {}
        for z, b in self._bins.items():
            model_params = self.parametrization.params_for_redshift(
                named_params, z
            )
            for d in b["datasets"]:
                model = d.predict(b["obs"], model_params)
                out[d.label] = {
                    "x": d.x,
                    "data": d.data,
                    "model": np.asarray(model),
                    "z": d.z,
                    "kind": d.kind,
                }
        return out

    def __getstate__(self):
        # Drop the live (unpicklable) Observables so each worker rebuilds them.
        state = self.__dict__.copy()
        state["_bins"] = None
        return state

    def __repr__(self):
        kinds = ", ".join(sorted({d.kind for d in self.datasets}))
        return (
            f"GaussianLikelihood({len(self.datasets)} datasets [{kinds}], "
            f"z={self.redshifts}, mode='{self.parametrization.mode}')"
        )


class _NonFinite(Exception):
    """Internal signal that a model produced a non-finite likelihood."""
