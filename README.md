# halogal

[![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python package for modeling UV luminosity functions and galaxy clustering using Halo Occupation Distribution (HOD) models for high-redshift galaxies.

## Features

- **Simple API**: Only redshift required, scientifically validated defaults
- **UV-Halo Mass Relation (UVHMR)**: Connect halo mass to UV luminosity through star formation
- **Halo Occupation Distribution**: Model central and satellite galaxy populations
- **Luminosity Functions**: Compute UV luminosity functions at high redshift
- **Galaxy Bias**: Calculate galaxy clustering bias
- **Correlation Functions**: Angular, real-space, and projected 2PCFs via halomod
- **Efficient Parameter Updates**: All observables accept inline `**params` for fast MCMC loops
- **Bayesian Inference**: `halogal.inference` fits the model to UVLF, clustering and number-density data with `emcee` or `dynesty` behind one interface
- **Dust Attenuation**: Self-consistent treatment following Bouwens+2013-14
- **Fitted Parameters**: Defaults from Shuntov+2025
- **Redshift Evolution**: Built-in parameter evolution with redshift

## Installation

### From PyPI (https://pypi.org/project/halogal/)

```bash
pip install halogal
```

### From source

```bash
git clone https://github.com/mshuntov/halogal.git
cd halogal
pip install -e .
```

To fit data you also need the samplers:

```bash
pip install "halogal[inference]"     # emcee, dynesty, corner, h5py, schwimmbad
```

The compiled dependencies (notably `camb`, pulled in by `halomod`) install most
reliably from conda-forge; `environment.yml` encodes a known-good stack:

```bash
micromamba env create -f environment.yml
micromamba activate halogal
pip install -e .
```

### Dependencies

- numpy >= 1.20, < 2 (required by `camb`)
- scipy >= 1.7
- astropy >= 5.0
- colossus >= 1.3
- halomod >= 1.4

## Quick Start

### Minimal Example

```python
import numpy as np
from halogal import HODModel
from halogal.model import Observables

# Create model — only redshift required!
# Uses fitted defaults from Shuntov+2025
model = HODModel(z=6.0)
obs = Observables(model)

# Compute luminosity function
MUV = np.linspace(-22, -16, 20)
phi = obs.luminosity_function(MUV)
```

### Galaxy Bias

```python
bias = obs.galaxy_bias(MUV)
```

### UV-Halo Mass Relation

```python
# UVHMR methods live directly on the model
Mh = 1e11  # M_sun
MUV = model.MUV(Mh)
sfr = model.sfr(Mh)

print(f"Halo mass {Mh:.2e} M_sun:")
print(f"  M_UV = {MUV:.2f}")
print(f"  SFR = {sfr:.1f} M_sun/yr")

# Inverse relation
Mh_recovered = model.Mhalo(MUV)
```

### Override Specific Parameters

```python
# At construction
model = HODModel(z=6.0, eps0=0.25, sigma_UV=0.5)
obs = Observables(model)

# Or inline when computing observables — no need to recreate anything
phi = obs.luminosity_function(MUV, eps0=0.3, sigma_UV=0.6)
```

### Efficient MCMC Fitting

All observable methods on `Observables` accept `**params` keyword arguments to
update the underlying model in-place before computing. This avoids object
re-creation in tight loops:

```python
obs = Observables(HODModel(z=6.0))
MUV = np.linspace(-22, -16, 20)

for eps0, sigma_UV in mcmc_samples:
    phi  = obs.luminosity_function(MUV, eps0=eps0, sigma_UV=sigma_UV)
    bg   = obs.galaxy_bias(MUV, eps0=eps0, sigma_UV=sigma_UV)
    mh   = obs.mean_halo_mass(-19.0, eps0=eps0, sigma_UV=sigma_UV)
    ngal = obs.number_density(-19.0, eps0=eps0, sigma_UV=sigma_UV)
```

For correlation functions, use the initialize/update pattern which leverages
halomod's internal caching:

```python
obs = Observables(HODModel(z=6.0))

# Initialize once (expensive)
result = obs.initialize_correlation_model(
    MUV_thresh1=-19.1, correlation_type='angular'
)

# Update efficiently in MCMC loop
for eps0, sigma_UV in mcmc_samples:
    result = obs.update_correlation_model(eps0=eps0, sigma_UV=sigma_UV)
    w_theta = result['correlation']
    theta   = result['separation']
```

### Fitting Data (`halogal.inference`)

```python
from halogal.inference import (
    Dataset, GaussianLikelihood, ParameterSet, Parameter,
    Uniform, ModelParametrization, EmceeSampler,
)

# What to fit: any mix of UVLF, clustering and number-density datasets
data = [
    Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=phi_err),
    Dataset.clustering(z=5.4, theta=theta_arcsec, w=w, sigma=w_err,
                       MUV_thresh1=-20.0,
                       corr_kwargs=dict(theta_min=0.5, theta_max=3600.0)),
]

# Which parameters are free (masses are sampled in log10)
params = ParameterSet([
    Parameter("eps0", Uniform(0.02, 0.4)),
    Parameter("logMc", Uniform(10.5, 12.5)),
    Parameter("sigma_UV", Uniform(0.2, 1.3)),
    Parameter("b", value=0.65, fixed=True),
])   # anything not listed stays at its published default

like = GaussianLikelihood(data, ModelParametrization("fixed_z", params))
like.check()   # evaluate once with all errors surfaced

result = EmceeSampler(params, like).run(
    nwalkers=32, processes=8, backend="chain.h5",
    convergence=True, max_steps=20000)
print(result)
result.corner_plot()
result.save("fit_z5.4")
```

Use `ModelParametrization("evolution", params)` with datasets at several
redshifts to fit the linear redshift evolution instead, and `DynestySampler`
for the Bayesian evidence. See `docs/inference.md`,
`examples/inference_tutorial.ipynb` and `examples/fit_shuntov25.py`.

### Redshift Evolution

```python
from halogal.models.parametrization import eps0_fz, Mc_fz
from halogal.config import DEFAULT_REDSHIFT_EVOLUTION

z_array = np.linspace(4, 8, 20)

# Get evolved parameters
eps0_z = eps0_fz(
    z_array, 
    deps_dz=DEFAULT_REDSHIFT_EVOLUTION['d_eps0_dz'],
    eps_off=DEFAULT_REDSHIFT_EVOLUTION['C_eps0']
)

Mc_z = 10**Mc_fz(
    z_array,
    dMc_dz=DEFAULT_REDSHIFT_EVOLUTION['d_logMc_dz'],
    Mc_off=DEFAULT_REDSHIFT_EVOLUTION['C_logMc']
)

# Compute with evolved parameters inline
obs = Observables(HODModel(z=z_array[0]))
for z, eps0, Mc in zip(z_array, eps0_z, Mc_z):
    phi = obs.luminosity_function(MUV, eps0=eps0, Mc=Mc)
```

### Compare to Observations

```python
from bouwens21_data import bouwens21, redshift_centers

data = bouwens21['z6']
z_obs = redshift_centers['z6']

obs = Observables(HODModel(z=z_obs))
MUV_model = np.linspace(-23, -15, 50)
phi_model = obs.luminosity_function(MUV_model)

plt.errorbar(data['M_AB'], data['Fi_k'], yerr=data['Fi_k_error'],
            fmt='o', label='Bouwens+2021')
plt.semilogy(MUV_model, phi_model, '-', label='Model')
plt.legend()
plt.show()
```

## Documentation

Full documentation is available at https://uvlf-hod.readthedocs.io/en/latest/#.

### Package Structure

```
halogal/
├── __init__.py          # Public API
├── config.py            # Configuration and defaults
├── cosmology.py         # Halo mass function and bias
├── model.py             # Unified UVHMR, HOD, and Observables
├── luminosity.py        # UV luminosity and dust
├── crosscorr.py         # Cross-correlations of two magnitude-selected samples
├── inference/           # MCMC / nested-sampling fits
│   ├── parameters.py       # Priors, Parameter, ParameterSet
│   ├── parametrization.py  # fixed-z / evolution parameter mapping
│   ├── data.py             # Dataset (UVLF, clustering, number density)
│   ├── likelihood.py       # GaussianLikelihood
│   ├── samplers.py         # EmceeSampler, DynestySampler
│   └── results.py          # InferenceResult
└── models/
    └── parametrization.py  # Redshift parametrizations
```

### Key Classes

- **`HODModel`**: Galaxy population model combining UVHMR + HOD (parameters and occupation functions)
- **`Observables`**: Compute observables (UVLF, bias, correlation functions) from an `HODModel`; supports efficient inline parameter updates via `**params`
- **`UVHMRModel`**: Base class for UV-halo mass relations only

## Model Architecture

```
UVHMRModel (base class)
├── Handles UV-halo mass relations
├── Methods: sfr(), MUV(), Mhalo(), star_formation_efficiency()
└── Parameters: z (required), eps0, Mc, a, b (optional)

HODModel (extends UVHMRModel)
├── Inherits all UVHMR methods
├── Adds occupation distributions: Ncen(), Nsat(), Ngal()
└── Additional parameters: sigma_UV, Mcut, Msat, asat (optional)

Observables (takes an HODModel)
├── luminosity_function(MUV, **params)
├── galaxy_bias(MUV, **params)
├── mean_halo_mass(MUV_thresh, **params)
├── mean_bias(MUV_thresh, **params)
├── number_density(MUV_thresh, **params)
├── initialize_correlation_model() / update_correlation_model()
└── compute_correlation_function()
```

All `Observables` methods accept `**params` (eps0, Mc, a, b, sigma_UV, Mcut, Msat, asat)
to update the model inline before computing.

## Default Parameters

All defaults from **Shuntov+2025** (2025A&A...699A.231S) at z~5.4:


### UVHMR Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `eps0` | 0.19 | Star formation efficiency |
| `Mc` | 10^11.64 M_☉ | Characteristic halo mass |
| `a` | 0.69 | Low-mass slope |
| `b` | 0.65 | High-mass slope |

### HOD Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sigma_UV` | 0.69 mag | UV magnitude scatter |
| `Mcut` | 10^9.57 M_☉ | Satellite cutoff mass |
| `Msat` | 10^12.65 M_☉ | Satellite normalization |
| `asat` | 0.85 | Satellite power-law slope |

### Redshift Evolution

All parameters evolve as: `param(z) = d_param/dz × z + C_param`

Evolution parameters available in `DEFAULT_REDSHIFT_EVOLUTION`.

## Examples

See the `examples/` directory for detailed examples:

- `basic_usage.ipynb`: Interactive Jupyter notebook with complete workflow
- `inference_tutorial.ipynb`: Injection-recovery tutorial for `halogal.inference`
- `fit_shuntov25.py`: Production-style parallel fit script (checkpointed)
- `bouwens21_data.py`: Observational data compilation

## Testing

Run tests with pytest:

```bash
pytest tests/ -m "not slow"   # fast unit tests (a few seconds)
pytest tests/                 # also runs the samplers (several minutes)
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new features
4. Submit a pull request

## Citation

If you use this package in your research, please cite:

```bibtex
@ARTICLE{2025A&A...699A.231S,
       author = {{Shuntov}, Marko and {Oesch}, Pascal A. and {Toft}, Sune and {Meyer}, Romain A. and {Covelo-Paz}, Alba and {Paquereau}, Louise and {Bouwens}, Rychard and {Brammer}, Gabriel and {Gelli}, Viola and {Giovinazzo}, Emma and {Herard-Demanche}, Thomas and {Illingworth}, Garth D. and {Mason}, Charlotte and {Naidu}, Rohan P. and {Weibel}, Andrea and {Xiao}, Mengyuan},
        title = "{Constraints on the early Universe star formation efficiency from galaxy clustering and halo modeling of H{\ensuremath{\alpha}} and [O III] emitters}",
      journal = {\aap},
     keywords = {galaxies: evolution, galaxies: high-redshift, galaxies: luminosity function, mass function, galaxies: statistics, Astrophysics of Galaxies},
         year = 2025,
        month = jul,
       volume = {699},
          eid = {A231},
        pages = {A231},
          doi = {10.1051/0004-6361/202554618},
archivePrefix = {arXiv},
       eprint = {2503.14280},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2025A&A...699A.231S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}


```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

Based on methodology from:
- Shuntov et al. 2025, A&A 699 A321
- Sabti et al. 2022
- Muñoz et al. 2023
- Bouwens et al. 2013, 2014
- And others...

## Contact

For questions or issues, please open an issue on GitHub or contact marko.shuntov@nbu.ku.dk.