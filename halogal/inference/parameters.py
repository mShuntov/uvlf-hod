"""Priors and parameter sets for MCMC inference.

This module decouples the flat *sampling vector* ``theta`` that a sampler
manipulates from the *named parameters* that the forward model consumes.

A :class:`Prior` knows how to evaluate its log-density (for ensemble samplers
such as emcee) **and** how to map a unit-cube coordinate to a parameter value
(the inverse-CDF / prior transform required by nested samplers such as dynesty).
Defining a prior once therefore drives both backends.

A :class:`ParameterSet` collects :class:`Parameter` specifications, tracks which
are free versus fixed, and provides the ``theta`` <-> named-dict mapping, the
total log-prior, the stacked unit-cube transform, and walker initialisation.
"""

import numpy as np
from scipy.special import erfinv

__all__ = [
    "Prior",
    "Uniform",
    "Normal",
    "LogUniform",
    "Parameter",
    "ParameterSet",
]


class Prior:
    """Abstract base class for a one-dimensional prior.

    Subclasses must implement :meth:`logpdf`, :meth:`unit_transform`,
    :meth:`ref`, and :meth:`sample`.
    """

    def logpdf(self, x):
        """Natural log of the prior probability density at ``x``."""
        raise NotImplementedError

    def unit_transform(self, u):
        """Map a unit-cube coordinate ``u`` in [0, 1] to a parameter value.

        This is the inverse cumulative distribution function, used as the
        prior transform for nested sampling.
        """
        raise NotImplementedError

    def ref(self):
        """Return a representative central value (for walker initialisation)."""
        raise NotImplementedError

    def sample(self, rng, size=None):
        """Draw random variate(s) from the prior."""
        raise NotImplementedError


class Uniform(Prior):
    """Uniform prior on the interval ``[low, high]``.

    Parameters
    ----------
    low, high : float
        Lower and upper bounds (``high`` must exceed ``low``).
    """

    def __init__(self, low, high):
        if not high > low:
            raise ValueError(f"Uniform requires high > low, got {low}, {high}")
        self.low = float(low)
        self.high = float(high)
        self._lognorm = -np.log(self.high - self.low)

    def logpdf(self, x):
        if self.low <= x <= self.high:
            return self._lognorm
        return -np.inf

    def unit_transform(self, u):
        return self.low + u * (self.high - self.low)

    def ref(self):
        return 0.5 * (self.low + self.high)

    def sample(self, rng, size=None):
        return rng.uniform(self.low, self.high, size=size)

    def __repr__(self):
        return f"Uniform(low={self.low}, high={self.high})"


class Normal(Prior):
    """Gaussian prior with mean ``mu`` and standard deviation ``sigma``.

    Parameters
    ----------
    mu : float
        Mean.
    sigma : float
        Standard deviation (must be positive).
    """

    def __init__(self, mu, sigma):
        if not sigma > 0:
            raise ValueError(f"Normal requires sigma > 0, got {sigma}")
        self.mu = float(mu)
        self.sigma = float(sigma)
        self._lognorm = -np.log(self.sigma * np.sqrt(2.0 * np.pi))

    def logpdf(self, x):
        return self._lognorm - 0.5 * ((x - self.mu) / self.sigma) ** 2

    def unit_transform(self, u):
        # Inverse CDF of the normal distribution via the error function.
        return self.mu + self.sigma * np.sqrt(2.0) * erfinv(2.0 * u - 1.0)

    def ref(self):
        return self.mu

    def sample(self, rng, size=None):
        return rng.normal(self.mu, self.sigma, size=size)

    def __repr__(self):
        return f"Normal(mu={self.mu}, sigma={self.sigma})"


class LogUniform(Prior):
    """Log-uniform (Jeffreys) prior on ``[low, high]`` with ``low > 0``.

    The density is proportional to ``1 / x``; equivalently the base-10
    logarithm of the parameter is uniform. Useful for a positive scale
    parameter sampled in linear space. (When a mass is sampled directly in
    ``log10`` space, a :class:`Uniform` prior on that log quantity is the
    natural choice instead.)

    Parameters
    ----------
    low, high : float
        Positive bounds with ``high`` > ``low``.
    """

    def __init__(self, low, high):
        if not (low > 0 and high > low):
            raise ValueError(
                f"LogUniform requires 0 < low < high, got {low}, {high}"
            )
        self.low = float(low)
        self.high = float(high)
        self._log_low = np.log10(self.low)
        self._log_high = np.log10(self.high)
        self._lognorm = -np.log(np.log(self.high / self.low))

    def logpdf(self, x):
        if self.low <= x <= self.high:
            return self._lognorm - np.log(x)
        return -np.inf

    def unit_transform(self, u):
        return 10.0 ** (self._log_low + u * (self._log_high - self._log_low))

    def ref(self):
        return 10.0 ** (0.5 * (self._log_low + self._log_high))

    def sample(self, rng, size=None):
        log_x = rng.uniform(self._log_low, self._log_high, size=size)
        return 10.0 ** log_x

    def __repr__(self):
        return f"LogUniform(low={self.low}, high={self.high})"


class Parameter:
    """Specification of a single model parameter.

    Parameters
    ----------
    name : str
        Parameter name. Must match the sampling-space names expected by the
        :class:`~halogal.inference.parametrization.ModelParametrization`
        (e.g. ``"eps0"``, ``"logMc"`` in fixed-z mode, or ``"d_eps0_dz"``,
        ``"C_eps0"`` in evolution mode).
    prior : Prior, optional
        Prior on the parameter. Required for free parameters; ignored for
        fixed ones.
    fixed : bool, optional
        If True the parameter is held at ``value`` and not sampled.
    value : float, optional
        For a fixed parameter, the held value (required). For a free
        parameter, an optional initial guess used to seed walkers; if omitted
        the prior's reference value is used.
    """

    def __init__(self, name, prior=None, fixed=False, value=None):
        self.name = name
        self.fixed = bool(fixed)
        self.value = value
        if self.fixed:
            if value is None:
                raise ValueError(f"Fixed parameter '{name}' requires a value")
            self.prior = None
        else:
            if prior is None:
                raise ValueError(f"Free parameter '{name}' requires a prior")
            self.prior = prior

    def initial_center(self):
        """Central value used to seed walkers (initial guess or prior ref)."""
        if self.value is not None:
            return float(self.value)
        return float(self.prior.ref())

    def __repr__(self):
        if self.fixed:
            return f"Parameter('{self.name}', fixed={self.value})"
        return f"Parameter('{self.name}', {self.prior!r})"


class ParameterSet:
    """An ordered collection of :class:`Parameter` objects.

    The free parameters define the sampling vector ``theta`` (in the order they
    were supplied). Fixed parameters never enter ``theta`` but are always
    present in the named dictionary returned by :meth:`to_dict`.

    Parameters
    ----------
    parameters : list of Parameter
        The parameter specifications.
    """

    def __init__(self, parameters):
        self.parameters = list(parameters)
        names = [p.name for p in self.parameters]
        if len(names) != len(set(names)):
            raise ValueError("Duplicate parameter names in ParameterSet")
        self.free = [p for p in self.parameters if not p.fixed]
        self.fixed_params = [p for p in self.parameters if p.fixed]
        self.free_names = [p.name for p in self.free]

    @property
    def ndim(self):
        """Number of free (sampled) parameters."""
        return len(self.free)

    @property
    def names(self):
        """Names of the free parameters, in sampling order."""
        return list(self.free_names)

    def to_dict(self, theta):
        """Map a free-parameter vector ``theta`` to a full named dict.

        Parameters
        ----------
        theta : array_like
            Values of the free parameters, in sampling order.

        Returns
        -------
        dict
            All parameters (free and fixed) keyed by name.
        """
        theta = np.asarray(theta)
        if theta.shape[-1] != self.ndim:
            raise ValueError(
                f"theta has length {theta.shape[-1]}, expected {self.ndim}"
            )
        params = {p.name: float(p.value) for p in self.fixed_params}
        for name, val in zip(self.free_names, theta):
            params[name] = float(val)
        return params

    def log_prior(self, theta):
        """Total log-prior of the free parameters at ``theta``."""
        lp = 0.0
        for p, val in zip(self.free, theta):
            lp += p.prior.logpdf(val)
            if not np.isfinite(lp):
                return -np.inf
        return lp

    def unit_transform(self, u):
        """Map a unit-cube vector ``u`` to ``theta`` (dynesty prior transform).

        Parameters
        ----------
        u : array_like
            Coordinates in the unit hypercube, one per free parameter.

        Returns
        -------
        ndarray
            Corresponding free-parameter values.
        """
        u = np.asarray(u)
        return np.array(
            [p.prior.unit_transform(ui) for p, ui in zip(self.free, u)]
        )

    def initial(self, nwalkers, rng=None, scatter=1e-3):
        """Initial walker positions for an ensemble sampler.

        A small Gaussian ball is placed around each parameter's central value
        (its initial guess, or the prior reference), and clipped to remain
        inside any bounded prior support so every walker starts with finite
        log-probability.

        Parameters
        ----------
        nwalkers : int
            Number of walkers.
        rng : numpy.random.Generator, optional
            Random generator. A default is created if omitted.
        scatter : float, optional
            Fractional scale of the initial ball. Default 1e-3.

        Returns
        -------
        ndarray, shape (nwalkers, ndim)
            Initial positions.
        """
        rng = np.random.default_rng() if rng is None else rng
        centers = np.array([p.initial_center() for p in self.free])
        scales = np.where(np.abs(centers) > 0, np.abs(centers) * scatter, scatter)
        pos = centers[None, :] + scales[None, :] * rng.standard_normal(
            (nwalkers, self.ndim)
        )
        # Clip into bounded prior support (open interval) to keep logp finite.
        for j, p in enumerate(self.free):
            prior = p.prior
            if isinstance(prior, Uniform):
                lo, hi = prior.low, prior.high
            elif isinstance(prior, LogUniform):
                lo, hi = prior.low, prior.high
            else:
                continue
            eps = 1e-9 * (hi - lo)
            pos[:, j] = np.clip(pos[:, j], lo + eps, hi - eps)
        return pos

    def __len__(self):
        return self.ndim

    def __repr__(self):
        return (
            f"ParameterSet({self.ndim} free, "
            f"{len(self.fixed_params)} fixed: {self.free_names})"
        )
