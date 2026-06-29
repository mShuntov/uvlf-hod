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


# --- Worker-global state for pool-based parallel emcee --------------------
# Building the per-redshift-bin Observables (the halomod correlation model) is
# expensive (~1 s). A naive multiprocessing pool re-pickles the likelihood on
# every emcee step, and because GaussianLikelihood.__getstate__ drops the live
# (unpicklable) halomod state, each worker would rebuild it on *every*
# evaluation -- erasing the parallel speedup. Instead we build it ONCE per
# worker via a Pool initializer and read it from this module global; the emcee
# step then only ships the cheap, module-level log-prob below.
_WORKER_STATE = {}


def _init_pool_worker(parameter_set, likelihood):
    """Pool initializer: store and build the likelihood once per worker."""
    _WORKER_STATE["parameter_set"] = parameter_set
    _WORKER_STATE["likelihood"] = likelihood
    likelihood._ensure_setup()  # build the halomod model once in this process


def _pool_worker_log_prob(theta):
    """Module-level log-posterior using the per-worker likelihood."""
    ps = _WORKER_STATE["parameter_set"]
    like = _WORKER_STATE["likelihood"]
    lp = ps.log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + like.log_likelihood(ps.to_dict(theta))


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

    def run(self, nwalkers=None, nsteps=2000, processes=None, pool=None,
            progress=True, burnin=None, thin=1, initial=None, seed=None,
            moves=None):
        """Run the ensemble sampler.

        Parameters
        ----------
        nwalkers : int, optional
            Number of walkers. Defaults to ``max(2 * ndim + 2, 32)``.
        nsteps : int, optional
            Number of steps per walker. Default 2000.
        processes : int, optional
            If given (>1), run in parallel across this many local processes.
            This is the recommended way to parallelise a clustering fit: the
            likelihood (including the expensive halomod model) is built **once
            per worker** via a pool initializer, so the ~Ncore speedup is
            actually realised. Mutually exclusive with ``pool``.
        pool : object, optional
            An explicit pool with a ``map`` method (e.g. a ``schwimmbad`` MPI
            pool for cluster runs). Note: with an external pool the likelihood
            is re-sent each step and rebuilt per evaluation, so prefer
            ``processes`` on a single machine. Mutually exclusive with
            ``processes``.
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

        if processes is not None and pool is not None:
            raise ValueError("Pass either 'processes' or 'pool', not both")

        ndim = self.parameter_set.ndim
        if nwalkers is None:
            nwalkers = max(2 * ndim + 2, 32)
        if burnin is None:
            burnin = nsteps // 3

        rng = np.random.default_rng(seed)
        if initial is None:
            initial = self.parameter_set.initial(nwalkers, rng=rng)

        # Choose the log-prob callable and (optionally) build an owned pool that
        # constructs the likelihood once per worker.
        own_pool = None
        if processes is not None and processes > 1:
            from multiprocessing import Pool
            own_pool = Pool(
                processes,
                initializer=_init_pool_worker,
                initargs=(self.parameter_set, self.likelihood),
            )
            pool = own_pool
            log_prob_fn = _pool_worker_log_prob
        else:
            log_prob_fn = self._log_prob

        try:
            sampler = emcee.EnsembleSampler(
                nwalkers, ndim, log_prob_fn, pool=pool, moves=moves
            )
            sampler.run_mcmc(initial, nsteps, progress=progress)
        finally:
            if own_pool is not None:
                own_pool.close()
                own_pool.join()

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
