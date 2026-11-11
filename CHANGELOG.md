# Changelog

## 0.2.0 (unreleased)

### Added

- **`halogal.inference`**: Bayesian fitting of the UVHMR/HOD model to UV
  luminosity functions, angular/projected clustering and number densities.
  - `Parameter`, `ParameterSet` and the `Uniform`, `Normal`, `LogUniform`
    priors; any parameter can be free or fixed.
  - `ModelParametrization` with two regimes: `fixed_z` (HOD parameters at one
    redshift) and `evolution` (linear redshift-evolution coefficients fitted
    jointly across bins).
  - `Dataset.uvlf`, `Dataset.clustering`, `Dataset.number_density`, each with a
    full covariance or diagonal errors.
  - `GaussianLikelihood`, which keeps one `halomod` model per redshift bin and
    only updates it while sampling; `check()` validates a setup before a run.
  - `EmceeSampler` and `DynestySampler` behind one interface. Both take
    `processes=N` and build the model once per worker.
  - `EmceeSampler.run`: HDF5 checkpoint/resume (`backend=`) and
    autocorrelation-based stopping (`convergence=True`).
  - `InferenceResult`: summaries, weighted quantiles, corner plot, evidence,
    `save()` / `load()`.
- `halogal.crosscorr`: `Selection` and `CrossCorrelation` for cross-correlations
  between two magnitude-selected samples.
- `environment.yml` with a known-good conda-forge stack, and the
  `halogal[inference]` extra.

### Fixed

- Cosmology is defined once (`config.DEFAULT_COSMO_PARAMS`) and passed to
  `halomod` as well as Colossus; previously the clustering ran with `halomod`'s
  own `sigma8`/`ns`/transfer function while the UVLF used halogal's.
- The `halomod` halo mass function is converted to `dn/dlog10M` like the
  Colossus one; absolute quantities on that path were low by `1/ln(10)`.
- `h` in the unit conversions follows the active cosmology instead of a
  hard-coded 0.7.
- Integer magnitude or mass arrays no longer truncate the UVLF / `N_sat`.
- Misspelled parameter names passed to `Observables` methods raise `TypeError`
  instead of being ignored.
- Magnitude-binned clustering samples (`MUV_thresh2`) are handled correctly
  and mis-ordered thresholds are rejected.
- Two typos in the Bouwens+2021 UVLF errors shipped in
  `examples/bouwens21_data.py` (z~5 and z~6, a factor of 10 too small).

### Changed

- Requires Python >= 3.9 and `numpy < 2` (a `camb` limitation).
- Only the `halogal` package is installed; the top-level `tests` and `examples`
  directories are no longer shipped as packages.
- `tests/test_core.py`, which targeted a module that no longer exists, is
  replaced by `tests/test_model.py`.

## 0.1.0

Initial release: UV-halo mass relation, HOD, UV luminosity function, galaxy
bias and `halomod` correlation functions.
