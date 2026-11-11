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
            moves=None, backend=None, convergence=False, max_steps=None,
            check_every=1000, n_tau=50, tau_change=0.01):
        """Run the ensemble sampler.

        Parameters
        ----------
        nwalkers : int, optional
            Number of walkers. Defaults to ``max(2 * ndim + 2, 32)``.
        nsteps : int, optional
            Number of steps per walker (fixed-length mode). Default 2000.
            Ignored when ``convergence=True`` (use ``max_steps`` instead).
        processes : int, optional
            If given (>1), run in parallel across this many local processes.
            The likelihood is built **once per worker** via a pool
            initializer, so the ~Ncore speedup is realised. Mutually
            exclusive with ``pool``.
        pool : object, optional
            An explicit pool with a ``map`` method (e.g. a ``schwimmbad`` MPI
            pool). With an external pool the likelihood is rebuilt per
            evaluation, so prefer ``processes`` on a single machine.
        progress : bool, optional
            Show a progress bar.
        burnin : int, optional
            Steps to discard when flattening. Defaults to ``nsteps // 3``
            (fixed mode) or ``int(2 * max(tau))`` when converged.
        thin : int, optional
            Thinning factor for the flattened chain. Default 1.
        initial : ndarray, optional
            Initial walker positions ``(nwalkers, ndim)``. Defaults to a small
            ball from :meth:`ParameterSet.initial`. Ignored when resuming a
            non-empty ``backend``.
        seed : int, optional
            Seed for walker initialisation.
        moves : optional
            ``emcee`` moves specification.
        backend : str or emcee.backends.HDFBackend, optional
            Enables crash-safe, resumable chains. Pass a filename (an
            ``HDFBackend`` is opened on it) or a backend instance. If the
            backend already holds iterations, the run **resumes** from the
            stored state (``initial`` is ignored) -- the basis for restarting
            a walltime-limited cluster job. A fresh backend is ``reset`` to
            ``(nwalkers, ndim)``.
        convergence : bool, optional
            If True, sample until the integrated autocorrelation time ``tau``
            is stable rather than for a fixed ``nsteps``: every
            ``check_every`` steps, stop when the chain is longer than
            ``n_tau * tau`` for every parameter **and** ``tau`` changed by
            less than ``tau_change`` since the last check (the canonical
            emcee criterion). Caps at ``max_steps``.
        max_steps : int, optional
            Iteration cap in convergence mode. Defaults to ``nsteps``.
        check_every : int, optional
            Steps between convergence checks. Default 1000.
        n_tau : int, optional
            Require ``chain length > n_tau * tau``. Default 50.
        tau_change : float, optional
            Require the fractional change in ``tau`` below this. Default 0.01.

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

        # Optional HDF backend (checkpoint / resume).
        if isinstance(backend, str):
            backend = emcee.backends.HDFBackend(backend)
        resuming = False
        if backend is not None:
            try:
                resuming = backend.iteration > 0
            except (AttributeError, OSError):
                resuming = False
            if resuming:
                nwalkers, ndim_b = backend.shape
                if ndim_b != ndim:
                    raise ValueError(
                        f"backend ndim {ndim_b} != parameter ndim {ndim}"
                    )
            else:
                backend.reset(nwalkers, ndim)

        rng = np.random.default_rng(seed)
        if resuming:
            initial = None  # emcee continues from the backend's last state
        elif initial is None:
            initial = self.parameter_set.initial(nwalkers, rng=rng)
        else:
            initial = self._check_initial(initial, nwalkers, ndim)

        # Build the model once in this process *before* any worker exists.
        # Setup errors (e.g. a bad halomod configuration) then surface here
        # with a clean traceback. Raised inside a Pool initializer instead,
        # they kill the worker, multiprocessing respawns it, and the run hangs
        # forever without ever reporting the error.
        self.likelihood._ensure_setup()

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

        converged = False
        finished = False
        try:
            sampler = emcee.EnsembleSampler(
                nwalkers, ndim, log_prob_fn, pool=pool, moves=moves,
                backend=backend,
            )
            if convergence:
                converged = self._run_to_convergence(
                    sampler, initial, max_steps or nsteps, check_every,
                    n_tau, tau_change, progress,
                )
            else:
                sampler.run_mcmc(initial, nsteps, progress=progress)
            finished = True
        finally:
            if own_pool is not None:
                # On an error / Ctrl-C do not wait for queued evaluations.
                if finished:
                    own_pool.close()
                else:
                    own_pool.terminate()
                own_pool.join()

        try:
            autocorr = sampler.get_autocorr_time(tol=0)
        except Exception:
            autocorr = None

        if burnin is None:
            if converged and autocorr is not None:
                burnin = int(2 * np.max(autocorr))
            else:
                burnin = sampler.iteration // 3
        burnin = min(burnin, max(sampler.iteration - 1, 0))

        flat = sampler.get_chain(discard=burnin, thin=thin, flat=True)
        log_prob = sampler.get_log_prob(discard=burnin, thin=thin, flat=True)

        extras = {
            "autocorr_time": autocorr,
            "acceptance_fraction": float(np.mean(sampler.acceptance_fraction)),
            "log_prob": log_prob,
            "burnin": burnin,
            "thin": thin,
            "converged": converged,
            "n_iterations": int(sampler.iteration),
            "sampler": sampler,
        }
        return InferenceResult(
            names=self.parameter_set.names,
            samples=flat,
            parameter_set=self.parameter_set,
            extras=extras,
            backend="emcee",
        )

    def _check_initial(self, initial, nwalkers, ndim):
        """Validate user-supplied walker positions (shape and prior support)."""
        initial = np.asarray(initial, dtype=float)
        if initial.shape != (nwalkers, ndim):
            raise ValueError(
                f"initial has shape {initial.shape}, expected "
                f"(nwalkers, ndim) = ({nwalkers}, {ndim})"
            )
        outside = [
            i for i, theta in enumerate(initial)
            if not np.isfinite(self.parameter_set.log_prior(theta))
        ]
        if outside:
            raise ValueError(
                f"{len(outside)} of {nwalkers} initial walker positions lie "
                f"outside the prior support (first: walker {outside[0]}). "
                "Such walkers start at log-probability -inf and never move."
            )
        return initial

    @staticmethod
    def _run_to_convergence(sampler, initial, max_steps, check_every,
                            n_tau, tau_change, progress):
        """Sample until tau is stable (canonical emcee criterion).

        Returns True if convergence was reached before ``max_steps``.
        """
        # When resuming from a backend, ``initial`` is None; the low-level
        # ``sample`` generator (unlike ``run_mcmc``) needs the explicit last
        # state rather than None.
        if initial is None:
            initial = sampler.get_last_sample()
        old_tau = np.inf
        for _ in sampler.sample(initial, iterations=max_steps,
                                progress=progress):
            if sampler.iteration % check_every:
                continue
            tau = sampler.get_autocorr_time(tol=0)
            if not np.all(np.isfinite(tau)):
                continue
            long_enough = np.all(n_tau * tau < sampler.iteration)
            stable = np.all(np.abs(old_tau - tau) / tau < tau_change)
            old_tau = tau
            if long_enough and stable:
                return True
        return False


class DynestySampler(_BaseSampler):
    """Static nested sampler (``dynesty``).

    Also gradient-free; robust to multimodality and -- unlike emcee -- returns
    the Bayesian evidence ``log Z`` for model comparison.
    """

    def run(self, nlive=500, processes=None, pool=None, queue_size=None,
            dlogz=0.01, maxiter=None, maxcall=None, print_progress=True,
            seed=None, **kwargs):
        """Run nested sampling to a stopping criterion.

        Parameters
        ----------
        nlive : int, optional
            Number of live points. Default 500.
        processes : int, optional
            If given (>1), run in parallel across this many local processes
            using ``dynesty.pool.Pool``, which ships the likelihood to each
            worker **once**; the expensive halomod model is then built once
            per worker rather than on every evaluation. Mutually exclusive
            with ``pool``.
        pool : object, optional
            An explicit pool (``multiprocessing`` or ``schwimmbad``). When
            given, ``queue_size`` should match the number of workers; it is
            inferred from ``pool._processes`` when possible. With an external
            pool the likelihood is re-sent and rebuilt per evaluation, so
            prefer ``processes`` on a single machine.
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

        if processes is not None and pool is not None:
            raise ValueError("Pass either 'processes' or 'pool', not both")

        ndim = self.parameter_set.ndim
        # Fail fast on setup errors, in this process (see EmceeSampler.run).
        self.likelihood._ensure_setup()

        loglike, ptform = self._log_likelihood, self._prior_transform
        own_pool = None
        if processes is not None and processes > 1:
            from dynesty.pool import Pool as DynestyPool
            own_pool = DynestyPool(processes, loglike, ptform)
            own_pool.__enter__()
            loglike, ptform = own_pool.loglike, own_pool.prior_transform
            pool, queue_size = own_pool, processes
        elif pool is not None and queue_size is None:
            queue_size = getattr(pool, "_processes", None) or getattr(
                pool, "size", None
            )

        rstate = np.random.default_rng(seed)
        finished = False
        try:
            sampler = dynesty.NestedSampler(
                loglike,
                ptform,
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
            finished = True
        finally:
            if own_pool is not None:
                own_pool.__exit__(None, None, None)  # terminate + join

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
