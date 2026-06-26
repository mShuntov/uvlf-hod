"""Unified container for posterior samples from either backend."""

import numpy as np

__all__ = ["InferenceResult"]


def _weighted_quantile(values, quantiles, weights=None):
    """Weighted quantiles of a 1-D array (unweighted if ``weights`` is None)."""
    values = np.asarray(values, dtype=float)
    quantiles = np.asarray(quantiles, dtype=float)
    if weights is None:
        return np.quantile(values, quantiles)
    weights = np.asarray(weights, dtype=float)
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cdf = np.cumsum(weights) - 0.5 * weights
    cdf /= np.sum(weights)
    return np.interp(quantiles, cdf, values)


class InferenceResult:
    """Posterior samples plus summary, evidence, and plotting helpers.

    Parameters
    ----------
    names : list of str
        Free-parameter names, in column order of ``samples``.
    samples : ndarray, shape (nsamples, ndim)
        Posterior samples.
    parameter_set : ParameterSet, optional
        The parameter set sampled; lets :meth:`median_params` reconstruct the
        full named dict (including fixed parameters).
    weights : ndarray, optional
        Per-sample weights (dynesty). Equal weights assumed if omitted.
    log_evidence, log_evidence_err : float, optional
        Bayesian log-evidence and its uncertainty (dynesty).
    extras : dict, optional
        Backend-specific diagnostics (e.g. ``autocorr_time``,
        ``acceptance_fraction``, the raw ``sampler``).
    backend : str, optional
        Name of the sampler that produced the result.
    """

    def __init__(self, names, samples, parameter_set=None, weights=None,
                 log_evidence=None, log_evidence_err=None, extras=None,
                 backend=None):
        self.names = list(names)
        self.samples = np.atleast_2d(np.asarray(samples, dtype=float))
        self.parameter_set = parameter_set
        self.weights = None if weights is None else np.asarray(weights, float)
        self.log_evidence = log_evidence
        self.log_evidence_err = log_evidence_err
        self.extras = extras or {}
        self.backend = backend

    def get_chain(self):
        """Return the flat sample array, shape (nsamples, ndim)."""
        return self.samples

    def summary(self, quantiles=(0.16, 0.50, 0.84)):
        """Per-parameter median and credible interval.

        Returns
        -------
        dict
            ``{name: {"median": m, "lower": m - q16, "upper": q84 - m,
            "q16": .., "q50": .., "q84": ..}}`` using the supplied quantiles.
        """
        ql, qm, qh = quantiles
        out = {}
        for j, name in enumerate(self.names):
            lo, med, hi = _weighted_quantile(
                self.samples[:, j], [ql, qm, qh], self.weights
            )
            out[name] = {
                "median": med,
                "lower": med - lo,
                "upper": hi - med,
                "q16": lo,
                "q50": med,
                "q84": hi,
            }
        return out

    def median_theta(self):
        """Median of the free parameters as a vector (weighted if applicable)."""
        return np.array(
            [
                _weighted_quantile(self.samples[:, j], [0.5], self.weights)[0]
                for j in range(len(self.names))
            ]
        )

    def median_params(self):
        """Full named parameter dict at the posterior median.

        Includes fixed parameters when a ``parameter_set`` was supplied, so the
        result can be passed straight to ``GaussianLikelihood.predict_all``.
        """
        theta = self.median_theta()
        if self.parameter_set is not None:
            return self.parameter_set.to_dict(theta)
        return dict(zip(self.names, theta))

    def corner_plot(self, truths=None, **kwargs):
        """Corner plot of the posterior (requires ``corner``).

        Uses the project's ``marko.mplstyle`` when available. Extra keyword
        arguments are forwarded to ``corner.corner``.
        """
        import corner
        import matplotlib.pyplot as plt

        opts = dict(labels=self.names, show_titles=True,
                    quantiles=[0.16, 0.5, 0.84])
        if self.weights is not None:
            opts["weights"] = self.weights
        if truths is not None:
            opts["truths"] = truths
        opts.update(kwargs)

        try:
            with plt.style.context("marko.mplstyle"):
                return corner.corner(self.samples, **opts)
        except (OSError, IOError):
            return corner.corner(self.samples, **opts)

    def __str__(self):
        lines = [f"InferenceResult (backend={self.backend}, "
                 f"nsamples={len(self.samples)})"]
        if self.log_evidence is not None:
            err = "" if self.log_evidence_err is None \
                else f" +/- {self.log_evidence_err:.3f}"
            lines.append(f"  log Z = {self.log_evidence:.3f}{err}")
        for name, s in self.summary().items():
            lines.append(
                f"  {name:>14s} = {s['median']:.4g} "
                f"(+{s['upper']:.3g} / -{s['lower']:.3g})"
            )
        return "\n".join(lines)

    __repr__ = __str__
