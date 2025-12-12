# uvlf-hod

[![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python package for modeling UV luminosity functions and galaxy clustering using Halo Occupation Distribution (HOD) models for high-redshift galaxies.

## Features

- **Simple API**: Only redshift required, scientifically validated defaults
- **UV-Halo Mass Relation (UVHMR)**: Connect halo mass to UV luminosity through star formation
- **Halo Occupation Distribution**: Model central and satellite galaxy populations
- **Luminosity Functions**: Compute UV luminosity functions at high redshift
- **Galaxy Bias**: Calculate galaxy clustering bias
- **Dust Attenuation**: Self-consistent treatment following Bouwens+2013-14
- **Fitted Parameters**: Defaults from Shuntov+2025
- **Redshift Evolution**: Built-in parameter evolution with redshift

## Installation

### From PyPI (when published)

```bash
pip install uvlf-hod
```

### From source

```bash
git clone https://github.com/mshuntov/uvlf-hod.git
cd uvlf-hod
pip install -e .
```

### Dependencies

- numpy >= 1.20
- scipy >= 1.7
- astropy >= 5.0
- colossus >= 1.3
- halomod >= 1.4

## Quick Start

### Minimal Example

```python
import numpy as np
import matplotlib.pyplot as plt
from uvlf_hod import HODModel

# Create model - only redshift required!
# Uses fitted defaults from Shuntov+2025
model = HODModel(z=6.0)

# Compute luminosity function
MUV = np.linspace(-22, -16, 20)
phi = model.luminosity_function(MUV)

# Plot
plt.semilogy(MUV, phi)
plt.xlabel('$M_{UV}$')
plt.ylabel('$\\Phi$ [Mpc$^{-3}$ mag$^{-1}$]')
plt.show()
```

### Galaxy Bias

```python
# Compute bias with the same model
bias = model.galaxy_bias(MUV)

plt.plot(MUV, bias)
plt.xlabel('$M_{UV}$')
plt.ylabel('Galaxy Bias $b_g$')
plt.show()
```

### UV-Halo Mass Relation

```python
# UVHMR methods are built-in
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
# Use defaults but change star formation efficiency
model = HODModel(z=6.0, eps0=0.25)

# Override multiple parameters
model = HODModel(z=6.0, eps0=0.2, Mc=10**12, sigma_UV=0.5)
```

### Redshift Evolution

```python
from uvlf_hod.models.parametrization import eps0_fz, Mc_fz
from uvlf_hod.config import DEFAULT_REDSHIFT_EVOLUTION

# Parameters evolve with redshift using fitted evolution
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

# Create models with evolved parameters
for z, eps0, Mc in zip(z_array, eps0_z, Mc_z):
    model = HODModel(z=z, eps0=eps0, Mc=Mc)
    phi = model.luminosity_function(MUV)
```

### Compare to Observations

```python
from bouwens21_data import bouwens21, redshift_centers

# Load data
obs = bouwens21['z6']
z_obs = redshift_centers['z6']

# Create model with defaults (fitted to this data!)
model = HODModel(z=z_obs)

# Compute and compare
MUV_model = np.linspace(-23, -15, 50)
phi_model = model.luminosity_function(MUV_model)

plt.errorbar(obs['M_AB'], obs['Fi_k'], yerr=obs['Fi_k_error'],
            fmt='o', label='Bouwens+2021')
plt.semilogy(MUV_model, phi_model, '-', label='Model')
plt.legend()
plt.show()
```

## Documentation

Full documentation is available at [link to docs].

### Package Structure

```
uvlf_hod/
├── __init__.py          # Public API
├── config.py            # Configuration and defaults
├── cosmology.py         # Halo mass function and bias
├── model.py             # Unified UVHMR and HOD models
├── luminosity.py        # UV luminosity and dust
└── models/
    └── parametrization.py  # Redshift parametrizations
```

### Key Classes

- **`HODModel`**: Main class for all galaxy population modeling (recommended)
- **`UVHMRModel`**: Base class for UV-halo mass relations only
- **`HaloMassFunction`**: Halo mass function with caching
- **`CosmologyConfig`**: Cosmology configuration

## Model Architecture

Clean, simple class hierarchy:

```
UVHMRModel (base class)
├── Handles UV-halo mass relations
├── Methods: sfr(), MUV(), Mhalo()
└── Parameters: z (required), eps0, Mc, a, b (optional)

HODModel (extends UVHMRModel)
├── Inherits all UVHMR methods
├── Adds occupation distributions
├── Methods: luminosity_function(), galaxy_bias(), Ncen(), Nsat()
└── Additional parameters: sigma_UV, Mcut, Msat, asat (optional)
```

## Default Parameters

<<<<<<< HEAD
All defaults from **Shuntov+2025** (2025A&A...699A.231S) z~5.4:
=======
All defaults from **Shuntov+2025** (2025A&A...699A.231S) at z~5.4:
>>>>>>> 2c62ade405dc58526b295d4dfca7fd65cd02f497

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
- `bouwens21_data.py`: Observational data compilation

## Testing

Run tests with pytest:

```bash
pytest tests/
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
