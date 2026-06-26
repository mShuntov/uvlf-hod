"""Observational data containers for inference.

A :class:`Dataset` bundles one observable measured at one redshift together
with its covariance, and knows how to obtain the corresponding model
prediction from a :class:`~halogal.model.Observables` instance. Three kinds are
supported:

``uvlf``
    UV luminosity function Phi(M_UV) at a grid of magnitudes.
``clustering``
    Angular (or projected/real) correlation function w(theta) of a
    magnitude-limited sample.
``number_density``
    Cumulative number density n(<M_UV thresh) (a scalar constraint).

The likelihood is Gaussian, so each dataset stores a covariance matrix. If only
a 1-D ``sigma`` is given it is promoted to a diagonal covariance -- the
covariance path is the same either way, so a full covariance (e.g. a jackknife
estimate for clustering) can be supplied later without code changes. The
Cholesky factor and Gaussian normalisation constant are cached at construction.
"""

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from .parametrization import MODEL_PARAM_NAMES

__all__ = ["Dataset"]

_KINDS = ("uvlf", "clustering", "number_density")


def _build_covariance(n, cov, sigma):
    """Return an (n, n) covariance from either ``cov`` or a 1-D ``sigma``."""
    if cov is not None and sigma is not None:
        raise ValueError("Provide either 'cov' or 'sigma', not both")
    if cov is None and sigma is None:
        raise ValueError("One of 'cov' or 'sigma' is required")
    if cov is not None:
        cov = np.atleast_2d(np.asarray(cov, dtype=float))
        if cov.shape != (n, n):
            raise ValueError(
                f"cov has shape {cov.shape}, expected ({n}, {n})"
            )
        return cov
    sigma = np.atleast_1d(np.asarray(sigma, dtype=float))
    if sigma.shape != (n,):
        raise ValueError(f"sigma has length {sigma.shape}, expected ({n},)")
    return np.diag(sigma ** 2)


class Dataset:
    """A single observable at a single redshift with its covariance.

    Use the classmethod constructors :meth:`uvlf`, :meth:`clustering`, and
    :meth:`number_density` rather than instantiating directly.
    """

    def __init__(self, kind, z, data, cov=None, sigma=None, label=None):
        if kind not in _KINDS:
            raise ValueError(f"kind must be one of {_KINDS}, got {kind!r}")
        self.kind = kind
        self.z = float(z)
        self.data = np.atleast_1d(np.asarray(data, dtype=float))
        self.label = label or f"{kind}_z{z:g}"

        self.cov = _build_covariance(self.data.size, cov, sigma)
        # Cache Cholesky factor and the Gaussian normalisation constant.
        self._cho = cho_factor(self.cov)
        logdet = 2.0 * np.sum(np.log(np.diag(self._cho[0])))
        self._lognorm = -0.5 * (self.data.size * np.log(2.0 * np.pi) + logdet)

        # kind-specific attributes (set by constructors)
        self.x = None
        self.MUV_thresh1 = None
        self.MUV_thresh2 = 0.0
        self.corr_kwargs = {}

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------
    @classmethod
    def uvlf(cls, z, MUV, phi, cov=None, sigma=None, label=None):
        """UV luminosity function data.

        Parameters
        ----------
        z : float
            Redshift.
        MUV : array_like
            Absolute UV magnitudes (data abscissa).
        phi : array_like
            Measured Phi(M_UV) in Mpc^-3 mag^-1.
        cov, sigma : array_like, optional
            Covariance (n, n) or 1-D Gaussian errors. Exactly one required.
        label : str, optional
            Human-readable label.
        """
        obj = cls("uvlf", z, phi, cov=cov, sigma=sigma, label=label)
        obj.x = np.atleast_1d(np.asarray(MUV, dtype=float))
        if obj.x.shape != obj.data.shape:
            raise ValueError("MUV and phi must have the same length")
        return obj

    @classmethod
    def number_density(cls, z, MUV_thresh, n, cov=None, sigma=None, label=None):
        """Cumulative number density n(<M_UV thresh).

        Parameters
        ----------
        z : float
            Redshift.
        MUV_thresh : float
            UV magnitude threshold (brighter than).
        n : float
            Measured number density in Mpc^-3.
        cov, sigma : optional
            1x1 covariance or scalar/1-element error.
        label : str, optional
            Label.
        """
        obj = cls("number_density", z, np.atleast_1d(n),
                  cov=cov, sigma=sigma, label=label)
        obj.MUV_thresh1 = float(MUV_thresh)
        return obj

    @classmethod
    def clustering(cls, z, theta, w, cov=None, sigma=None, MUV_thresh1=None,
                   MUV_thresh2=0.0, corr_kwargs=None, label=None):
        """Angular correlation function w(theta) of a magnitude-limited sample.

        Parameters
        ----------
        z : float
            Redshift.
        theta : array_like
            Angular separations in arcsec (data abscissa).
        w : array_like
            Measured w(theta).
        cov, sigma : array_like, optional
            Covariance (n, n) or 1-D errors. Exactly one required. For
            clustering a full covariance is strongly preferred because theta
            bins are correlated; a diagonal sigma is accepted but will
            underestimate parameter uncertainties.
        MUV_thresh1 : float
            UV magnitude threshold defining the clustered sample (required).
        MUV_thresh2 : float, optional
            Faint-end threshold for cross-correlations. Default 0.0.
        corr_kwargs : dict, optional
            Extra keyword arguments forwarded to
            :meth:`~halogal.model.Observables.initialize_correlation_model`
            (e.g. ``p1`` redshift distribution, ``correlation_type``,
            ``theta_min``/``theta_max``, ``znum``).
        label : str, optional
            Label.
        """
        if MUV_thresh1 is None:
            raise ValueError("clustering dataset requires MUV_thresh1")
        obj = cls("clustering", z, w, cov=cov, sigma=sigma, label=label)
        obj.x = np.atleast_1d(np.asarray(theta, dtype=float))
        if obj.x.shape != obj.data.shape:
            raise ValueError("theta and w must have the same length")
        obj.MUV_thresh1 = float(MUV_thresh1)
        obj.MUV_thresh2 = float(MUV_thresh2)
        obj.corr_kwargs = dict(corr_kwargs or {})
        return obj

    # ------------------------------------------------------------------
    # Prediction
    # ------------------------------------------------------------------
    def predict(self, observables, model_params):
        """Model prediction aligned with this dataset's abscissa.

        Parameters
        ----------
        observables : halogal.model.Observables
            Observables calculator for this dataset's redshift bin. For a
            clustering dataset its correlation model must already be initialised
            (the likelihood handles this).
        model_params : dict
            Model HOD parameters (linear masses), e.g. from
            :meth:`ModelParametrization.params_for_redshift`.

        Returns
        -------
        ndarray
            Model vector, same length as ``self.data``.
        """
        if self.kind == "uvlf":
            return np.atleast_1d(
                observables.luminosity_function(self.x, **model_params)
            )
        if self.kind == "number_density":
            n = observables.number_density(self.MUV_thresh1, **model_params)
            return np.atleast_1d(n)
        if self.kind == "clustering":
            corr_params = {
                k: v for k, v in model_params.items()
                if k in MODEL_PARAM_NAMES
            }
            result = observables.update_correlation_model(**corr_params)
            theta_model = np.asarray(result["separation"], dtype=float)
            w_model = np.asarray(result["correlation"], dtype=float)
            # Interpolate in log-separation (theta spans decades); keep w
            # linear so any sign is handled.
            return np.interp(
                np.log10(self.x), np.log10(theta_model), w_model
            )
        raise RuntimeError(f"Unknown kind {self.kind!r}")  # pragma: no cover

    def chi2(self, model):
        """Generalised chi-square ``r^T C^-1 r`` for a model vector."""
        r = self.data - np.asarray(model, dtype=float)
        return float(r @ cho_solve(self._cho, r))

    def log_likelihood(self, model):
        """Gaussian log-likelihood (with normalisation) for a model vector."""
        return self._lognorm - 0.5 * self.chi2(model)

    def __repr__(self):
        return f"Dataset(kind='{self.kind}', z={self.z:g}, n={self.data.size})"
