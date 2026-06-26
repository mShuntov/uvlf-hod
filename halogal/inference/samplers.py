"""Sampler backends: emcee (ensemble MCMC) and dynesty (nested sampling).

Both consume the same :class:`~halogal.inference.parameters.ParameterSet` and
:class:`~halogal.inference.likelihood.GaussianLikelihood` and return a common
:class:`~halogal.inference.results.InferenceResult`. The prior drives both: its
log-density feeds emcee's posterior, while its unit-cube transform serves as
dynesty's prior transform.

All callables passed to the samplers are bound methods (not closures) so that
they remain picklable when a ``pool`` is supplied for parallel evaluation.
"""

import numpy as np

from .results import InferenceResult

__all__ = ["EmceeSampler", "DynestySampler"]


class _BaseSampler:
    """Shared wiring between a parameter set and a likelihood."""

    def __init__(self, parameter_set, likelihood):
        self.parameter_set = parameter_set
        self.likelihood = likelihood
        # Ensure the parametrization's required names are all present.
        likelihood.parametrization.validate(parameter_set)

    def _log_likelihood(self, theta):
        named = self.parameter_set.to_dict(theta)
        return self.likelihood.log_likelihood(named)

    def _log_prob(self, theta):
        """Log-posterior up to a constant (for emcee)."""
        lp = self.parameter_set.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self._log_likelihood(theta)

    def _prior_transform(self, u):
        """Unit-cube -> parameter transform (for dynesty)."""
        return self.parameter_set.unit_transform(u)


class EmceeSampler(_BaseSampler):
    """Affine-invariant ensemble sampler (``emcee``).

    Gradient-free and embarrassingly parallel across walkers -- the right tool
    for this non-differentiable forward model, especially when clustering makes
    each likelihood call expensive.
    """

    def run(self, nwalkers=None, nsteps=2000, pool=None, progress=True,
            burnin=None, thin=1, initial=None, seed=None, moves=None):
        """Run the ensemble sampler.

        Parameters
        ----------
        nwalkers : int, optional
            Number of walkers. Defaults to ``max(2 * ndim + 2, 32)``.
        nsteps : int, optional
            Number of steps per walker. Default 2000.
        pool : object, optional
            Any object with a ``map`` method (``multiprocessing.Pool``,
            ``schwimmbad`` MPIPool, etc.) for parallel likelihood evaluation.
        progress : bool, optional
            Show a progress bar.
        burnin : int, optional
            Steps to discard. Defaults to ``nsteps // 3``.
        thin : int, optional
            Thinning factor for the flattened chain. Default 1.
        initial : ndarray, optional
            Initial walker positions, shape ``(nwalkers, ndim)``. Defaults to a
            small ball from :meth:`ParameterSet.initial`.
        seed : int, optional
            Seed for walker initialisation.
        moves : optional
            ``emcee`` moves specification.

        Returns
        -------
        InferenceResult
        """
        import emcee

        ndim = self.parameter_set.ndim
        if nwalkers is None:
            nwalkers = max(2 * ndim + 2, 32)
        if burnin is None:
            burnin = nsteps // 3

        rng = np.random.default_rng(seed)
        if initial is None:
            initial = self.parameter_set.initial(nwalkers, rng=rng)

        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, self._log_prob, pool=pool, moves=moves
        )
        sampler.run_mcmc(initial, nsteps, progress=progress)

        flat = sampler.get_chain(discard=burnin, thin=thin, flat=True)
        log_prob = sampler.get_log_prob(discard=burnin, thin=thin, flat=True)

        try:
            autocorr = sampler.get_autocorr_time(tol=0)
        except Exception:
            autocorr = None

        extras = {
            "autocorr_time": autocorr,
            "acceptance_fraction": float(np.mean(sampler.acceptance_fraction)),
            "log_prob": log_prob,
            "burnin": burnin,
            "thin": thin,
            "sampler": sampler,
        }
        return InferenceResult(
            names=self.parameter_set.names,
            samples=flat,
            parameter_set=self.parameter_set,
            extras=extras,
            backend="emcee",
        )


class DynestySampler(_BaseSampler):
    """Static nested sampler (``dynesty``).

    Also gradient-free; robust to multimodality and -- unlike emcee -- returns
    the Bayesian evidence ``log Z`` for model comparison.
    """

    def run(self, nlive=500, pool=None, queue_size=None, dlogz=0.01,
            maxiter=None, maxcall=None, print_progress=True, seed=None,
            **kwargs):
        """Run nested sampling to a stopping criterion.

        Parameters
        ----------
        nlive : int, optional
            Number of live points. Default 500.
        pool : object, optional
            Parallel pool (``multiprocessing`` or ``schwimmbad``). When given,
            ``queue_size`` should match the number of workers; it is inferred
            from ``pool._processes`` when possible.
        queue_size : int, optional
            Number of parallel likelihood evaluations.
        dlogz : float, optional
            Target evidence tolerance for stopping. Default 0.01.
        maxiter, maxcall : int, optional
            Optional hard caps on iterations / likelihood calls.
        print_progress : bool, optional
            Print the live progress summary.
        seed : int, optional
            Random seed.
        **kwargs
            Forwarded to ``dynesty.NestedSampler``.

        Returns
        -------
        InferenceResult
            With ``log_evidence`` / ``log_evidence_err`` populated and
            importance ``weights`` attached.
        """
        import dynesty

        ndim = self.parameter_set.ndim
        if pool is not None and queue_size is None:
            queue_size = getattr(pool, "_processes", None) or getattr(
                pool, "size", None
            )

        rstate = np.random.default_rng(seed)
        sampler = dynesty.NestedSampler(
            self._log_likelihood,
            self._prior_transform,
            ndim,
            nlive=nlive,
            pool=pool,
            queue_size=queue_size,
            rstate=rstate,
            **kwargs,
        )
        sampler.run_nested(
            dlogz=dlogz,
            maxiter=maxiter,
            maxcall=maxcall,
            print_progress=print_progress,
        )

        res = sampler.results
        logz = res.logz[-1]
        logzerr = res.logzerr[-1]
        # Normalised importance weights.
        weights = np.exp(res.logwt - logz)

        extras = {"results": res, "sampler": sampler}
        return InferenceResult(
            names=self.parameter_set.names,
            samples=res.samples,
            parameter_set=self.parameter_set,
            weights=weights,
            log_evidence=float(logz),
            log_evidence_err=float(logzerr),
            extras=extras,
            backend="dynesty",
        )
