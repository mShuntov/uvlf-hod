# `halogal.inference`: MCMC fitting of the UVHMR/HOD model

This document explains how the `halogal.inference` subpackage is built and how
the pieces fit together. The audience is the package author, so it assumes
familiarity with the forward model (`halogal.model.Observables`) and focuses on
the inference machinery layered on top of it.

All paths below are relative to `halogal/inference/` unless stated otherwise.

---

## 1. Purpose and design philosophy

The goal is to fit the UVHMR/HOD model to data (UV luminosity functions,
angular clustering, and number densities) and recover posteriors on the HOD
parameters.

The forward model is **not differentiable**: a single prediction runs NumPy
integrations, `scipy.special.erf` for the magnitude scatter, `colossus` for the
halo mass function, and `halomod` for the correlation function. There are no
gradients to exploit, so the inference uses **gradient-free samplers**: `emcee`
(affine-invariant ensemble MCMC) and `dynesty` (static nested sampling). Both
are exposed behind one interface.

The central design idea is to **decouple the flat sampling vector `theta`** that
a sampler manipulates **from the named parameters the model consumes**. A
sampler only ever sees an array of free parameters. A mapping layer turns that
array into a named dict, and a parametrization layer turns the named dict into
the per-redshift HOD parameters `Observables` expects. This decoupling is what
lets the *same* likelihood serve both a single fixed-redshift fit and a joint
fit of the redshift evolution, and lets the user freely choose which parameters
are free versus fixed without touching the model.

A secondary principle: **everything expensive is built once and cached**. The
Cholesky factor of each covariance is computed at `Dataset` construction; the
`halomod` correlation model is initialised once per redshift bin and only
*updated* inside the sampling loop.

---

## 2. Architecture: a pipeline of layers

```
   sampler (theta array)
        |
        v
  ParameterSet              parameters.py
   - log_prior(theta)        -> drives emcee
   - unit_transform(u)       -> drives dynesty
   - to_dict(theta)          -> named-parameter dict
        |
        v
  ModelParametrization      parametrization.py
   params_for_redshift(named, z)
   -> {eps0, Mc, a, b, sigma_UV, Mcut, Msat, asat}  (linear masses)
        |
        v
  GaussianLikelihood        likelihood.py
   - one Observables per z bin (cached halomod)
   - loops bins x datasets
        |
        v
  Dataset.predict(obs, params)   data.py
   -> model vector aligned to the data abscissa
        |
        v
  Dataset.chi2 / log_likelihood  (cached Cholesky)
        |
        v
  InferenceResult           results.py
   summary(), median_params(), corner_plot(), log_evidence
```

Each layer has one responsibility:

- **`parameters.py`** — `Prior` subclasses and the `ParameterSet`. Owns the
  `theta` <-> named-dict mapping, the total log-prior, the unit-cube transform,
  and walker initialisation.
- **`parametrization.py`** — `ModelParametrization`. Owns the choice of redshift
  regime and the `log10` -> linear mass conversion. The boundary between
  "sampling space" and "model space".
- **`data.py`** — `Dataset`. Owns one observable at one redshift, its
  covariance, the cached Cholesky factor, the forward `predict`, and the
  Gaussian `chi2` / `log_likelihood`.
- **`likelihood.py`** — `GaussianLikelihood`. Composes datasets into one total
  likelihood, holds the per-bin `Observables`, and handles the
  initialise/update caching for clustering.
- **`samplers.py`** — `EmceeSampler` and `DynestySampler`. Thin wiring that
  hands the right callables to each backend.
- **`results.py`** — `InferenceResult`. A unified container for posterior
  samples (weighted or not) plus summaries, evidence, and a corner plot.

The public API is re-exported from `halogal/inference/__init__.py`.

---

## 3. The two redshift modes

`ModelParametrization(mode, ...)` accepts `mode="fixed_z"` or
`mode="evolution"`. The canonical description of the eight HOD parameters lives
in one table, `HOD_PARAM_SPEC`:

```python
# (model_name, sample_name, is_log, slope_name, offset_name)
HOD_PARAM_SPEC = [
    ("eps0",     "eps0",     False, "d_eps0_dz",     "C_eps0"),
    ("Mc",       "logMc",    True,  "d_logMc_dz",    "C_logMc"),
    ("a",        "a",        False, "d_a_dz",        "C_a"),
    ("b",        "b",        False, "d_b_dz",        "C_b"),
    ("sigma_UV", "sigma_UV", False, "d_sigmaUV_dz",  "C_sigmaUV"),
    ("Mcut",     "logMcut",  True,  "d_Mcut_dz",     "C_Mcut"),
    ("Msat",     "logMsat",  True,  "d_Msat_dz",     "C_Msat"),
    ("asat",     "asat",     False, "d_asat_dz",     "C_asat"),
]
```

- `model_name` is the key `Observables` consumes (linear mass for the three
  mass parameters).
- `sample_name` is the fixed-z sampling-space name. Masses are sampled in
  `log10` (`logMc`, `logMcut`, `logMsat`).
- `is_log` marks the three mass parameters.
- `slope_name` / `offset_name` are the evolution-coefficient names, chosen to
  match the keys in `halogal.config.DEFAULT_REDSHIFT_EVOLUTION` so a fit can be
  compared to the published values directly.

`fixed_z_param_names()` returns the eight `sample_name`s; `evolution_param_names()`
returns the sixteen slope/offset coefficient names.

### `fixed_z` mode

The sampled quantities *are* the HOD parameters at a single redshift. For each
spec entry the code looks up `sample_name` and (for masses) raises `10**val`:

```python
val = lookup(sample_name)
if is_log:
    val = 10.0 ** val
out[model_name] = val
```

### `evolution` mode

The sampled quantities are the linear evolution coefficients. Per-bin values are
evaluated with `halogal.models.parametrization.RedshiftParametrization` — the
same engine behind the package's `*_fz` helper functions. With `use_linear=True`
(default) this is `param(z) = slope * z + offset`; with `use_linear=False` it is
the `z/(1+z)` normalised form.

```python
slope = lookup(slope_name)
offset = lookup(offset_name)
val = float(self._evolution.evaluate(z, slope, offset)[0])
if is_log:                 # mass coefficients act in log10 space
    val = 10.0 ** val
out[model_name] = val
```

### Conventions worth noting

- **The `log10` -> linear mass conversion lives only here**, at the boundary
  just before the model is called. Sampling space carries log masses everywhere;
  model space carries linear masses in M_sun.
- **Fallback to config defaults.** `lookup(name)` returns the value from the
  parameter dict if present, otherwise `self.defaults[name]`. The defaults come
  from `default_fixed_z_values()` (built from `DEFAULT_HOD_PARAMS`) or
  `default_evolution_values()` (a copy of `DEFAULT_REDSHIFT_EVOLUTION`). This is
  why a user can vary only a subset of parameters and leave the rest at their
  fiducial published values.
- **Validation.** `ModelParametrization.validate(parameter_set)` rejects names
  that are not recognised for the mode (catches typos) and confirms every
  required name is resolvable from either the parameter set or `defaults`.
  Samplers call this in `_BaseSampler.__init__`, so a bad name fails fast.

---

## 4. A single log-likelihood evaluation

Starting from a `theta` array (one entry per free parameter), the flow is:

1. **`theta -> named dict`.** `ParameterSet.to_dict(theta)` zips the free names
   onto `theta` and merges in the fixed parameters:

   ```python
   params = {p.name: float(p.value) for p in self.fixed_params}
   for name, val in zip(self.free_names, theta):
       params[name] = float(val)
   ```

2. **Per-bin model params.** `GaussianLikelihood._accumulate` loops over each
   redshift bin and calls
   `ModelParametrization.params_for_redshift(named, z)` to get the eight
   linear-mass HOD parameters for that `z`.

3. **Prediction.** For each dataset in the bin, `Dataset.predict(obs, params)`
   produces a model vector aligned to the dataset's abscissa:
   - `uvlf` -> `observables.luminosity_function(self.x, **params)`
   - `number_density` -> `observables.number_density(self.MUV_thresh1, **params)`
   - `clustering` -> `observables.update_correlation_model(**params)`, then
     interpolate `w(theta)` onto the data's `theta` grid in `log10(theta)`
     (linear in `w` so any sign is handled):

     ```python
     return np.interp(np.log10(self.x),
                      np.log10(theta_model), w_model)
     ```

4. **Gaussian chi-square with cached Cholesky.** Each `Dataset` stored
   `self._cho = cho_factor(self.cov)` and a Gaussian normalisation `_lognorm` at
   construction, so per-call cost is one back-substitution:

   ```python
   def chi2(self, model):
       r = self.data - model
       return float(r @ cho_solve(self._cho, r))

   def log_likelihood(self, model):
       return self._lognorm - 0.5 * self.chi2(model)
   ```

5. **Reduction.** `GaussianLikelihood.log_likelihood` sums the per-dataset
   log-likelihoods. If any prediction is non-finite it short-circuits to `-inf`
   (via an internal `_NonFinite` signal); if a model evaluation raises and
   `catch_errors=True` (the default), it also returns `-inf` rather than
   crashing the chain.

### Per-bin persistent `Observables` and clustering caching

`GaussianLikelihood._ensure_setup()` builds **one `Observables(HODModel(z))` per
unique redshift** and stores them in an `OrderedDict` keyed by `z`. For any bin
that has a clustering dataset it calls `initialize_correlation_model(...)`
**once** (the expensive `halomod` setup), forwarding `MUV_thresh1`,
`MUV_thresh2`, and the dataset's `corr_kwargs`. During sampling only
`update_correlation_model(**params)` runs, which reuses the cached `halomod`
model. This is the initialise-once / update-many pattern documented on
`Observables` itself.

Because the cache is keyed by a single `halomod` model per bin, **at most one
clustering dataset is allowed per redshift** — `GaussianLikelihood._validate`
enforces this.

---

## 5. Priors drive both samplers

A `Prior` defines its behaviour for *both* backends, so it is specified once:

- `logpdf(x)` — natural log prior density; feeds the emcee posterior.
- `unit_transform(u)` — inverse CDF mapping a unit-cube coordinate to a value;
  serves as the dynesty prior transform.
- `ref()` — a central value for walker initialisation.
- `sample(rng, size)` — random draws.

Three concrete priors are provided: `Uniform`, `Normal` (inverse CDF via
`erfinv`), and `LogUniform` (Jeffreys / `1/x`).

`ParameterSet` aggregates these across the free parameters:

```python
def log_prior(self, theta):          # emcee
    lp = 0.0
    for p, val in zip(self.free, theta):
        lp += p.prior.logpdf(val)
        if not np.isfinite(lp):
            return -np.inf
    return lp

def unit_transform(self, u):         # dynesty
    return np.array([p.prior.unit_transform(ui)
                     for p, ui in zip(self.free, u)])
```

In `samplers.py`, `_BaseSampler` wires these in: `_log_prob` (emcee) adds
`log_prior` to the likelihood and rejects out-of-support points early, while
`_prior_transform` (dynesty) is just `ParameterSet.unit_transform`.

`ParameterSet.initial(nwalkers)` seeds emcee with a small Gaussian ball around
each parameter's `initial_center()` (its `value` if given, else the prior
`ref()`), clipped just inside any bounded `Uniform`/`LogUniform` support so every
walker starts with finite log-probability.

---

## 6. Parallelism and the `__getstate__` rebuild hook

Both samplers accept a `pool` (anything with a `map`, e.g.
`multiprocessing.Pool` or a `schwimmbad` MPI pool). The callables handed to the
samplers are **bound methods** (`self._log_prob`, `self._log_likelihood`,
`self._prior_transform`), not closures, so they pickle cleanly to workers.

The catch is that the live per-bin `Observables` (and their `halomod` state) are
**not safely picklable**. `GaussianLikelihood` solves this by dropping them from
its pickled state:

```python
def __getstate__(self):
    state = self.__dict__.copy()
    state["_bins"] = None     # rebuilt lazily in each worker
    return state
```

When a worker first evaluates the likelihood, `_ensure_setup()` sees
`self._bins is None` and rebuilds the `Observables` (and re-initialises the
clustering model) in that process. `_ensure_setup` is idempotent, so this
happens once per worker.

`DynestySampler.run` additionally infers `queue_size` from the pool
(`pool._processes` or `pool.size`) when not given explicitly.

---

## 7. Minimal end-to-end example

A fixed-redshift fit combining a UVLF and a clustering measurement at z = 6:

```python
import numpy as np
from halogal.inference import (
    Dataset, GaussianLikelihood, ParameterSet, Parameter,
    Uniform, Normal, ModelParametrization, EmceeSampler,
)

# 1. Data: any mix of observables at this redshift.
data = [
    Dataset.uvlf(z=6, MUV=MUV, phi=phi, sigma=phi_err),
    Dataset.clustering(z=6, theta=th, w=w, cov=w_cov, MUV_thresh1=-19.1),
]

# 2. Free / fixed parameters (sampling-space names; masses in log10).
params = ParameterSet([
    Parameter("eps0",     Uniform(0.01, 0.5)),
    Parameter("logMc",    Normal(11.6, 0.5)),
    Parameter("a",        Uniform(0.3, 1.2)),
    Parameter("b",        value=0.65, fixed=True),
    Parameter("sigma_UV", Uniform(0.1, 1.5)),
    Parameter("logMcut",  Uniform(7.0, 12.0)),
    Parameter("logMsat",  Uniform(11.0, 14.0)),
    Parameter("asat",     Uniform(0.2, 1.5)),
])

# 3. Likelihood = data + parametrization.
like = GaussianLikelihood(data, ModelParametrization("fixed_z", params))

# 4. Sample, then read the posterior.
result = EmceeSampler(params, like).run(nwalkers=64, nsteps=5000)
print(result.summary())          # {name: {median, lower, upper, q16, q50, q84}}

# Overlay the best-fit model on the data:
best = result.median_params()    # full named dict incl. fixed params
predictions = like.predict_all(best)
```

For an evolution fit, supply `Parameter`s named after the evolution coefficients
(`"d_eps0_dz"`, `"C_eps0"`, ...), use `ModelParametrization("evolution", params)`,
and include datasets at several redshifts. To get the Bayesian evidence for
model comparison, swap in `DynestySampler(params, like).run(nlive=500)`; the
returned `InferenceResult` carries `log_evidence`, `log_evidence_err`, and
importance `weights` (which `summary` and `corner_plot` use automatically).

---

## 8. How to extend

**Add a new prior.** Subclass `Prior` and implement `logpdf`,
`unit_transform`, `ref`, and `sample`. That single class then works with both
backends. If it has a bounded support that walker initialisation should clip to,
add it to the `isinstance` check in `ParameterSet.initial`.

**Add a new observable kind.** Three coordinated changes in `data.py`:
1. Add the kind name to `_KINDS`.
2. Add a classmethod constructor (mirror `uvlf` / `clustering`) that stores the
   abscissa and any kind-specific attributes, and calls the base `__init__` with
   the measured vector plus `cov`/`sigma`.
3. Add a branch in `Dataset.predict` that calls the appropriate `Observables`
   method and returns a vector aligned to the dataset abscissa.

If the new kind needs an expensive one-time setup (as clustering does), also
extend `GaussianLikelihood._ensure_setup` to initialise it per bin, and consider
whether a per-bin uniqueness constraint belongs in `_validate`.

**Add a new sampler backend.** Subclass `_BaseSampler` in `samplers.py`. You
already inherit `_log_likelihood`, `_log_prob`, and `_prior_transform`. Implement
a `run(...)` method that drives your backend with the appropriate callable
(`_log_prob` for posterior-based samplers, `_log_likelihood` +
`_prior_transform` for likelihood/prior-transform samplers) and returns an
`InferenceResult` with `names`, `samples`, the `parameter_set`, and any backend
diagnostics in `extras` (plus `weights` / `log_evidence` if applicable). Keep
every sampler-facing callable a bound method so it stays picklable under a pool.
