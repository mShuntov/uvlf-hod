# uvlf-hod

[![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python package for modeling UV luminosity functions and galaxy clustering using Halo Occupation Distribution (HOD) models for high-redshift galaxies.

## Features

- **Unified Model Architecture**: Single class interface for all calculations
- **UV-Halo Mass Relation (UVHMR)**: Connect halo mass to UV luminosity through star formation
- **Halo Occupation Distribution**: Model central and satellite galaxy populations
- **Luminosity Functions**: Compute UV luminosity functions at high redshift
- **Galaxy Bias**: Calculate galaxy clustering bias
- **Dust Attenuation**: Self-consistent treatment following Bouwens+2013-14
- **Flexible Parametrization**: Redshift-dependent model parameters
- **Modern Python**: Clean API with class inheritance and comprehensive documentation

## Installation

### From PyPI (when published)

```bash
pip install uvlf-hod
```

### From source

```bash
git clone https://github.com/yourusername/uvlf-hod.git
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

### Basic UV Luminosity Function

```python
import numpy as np
import matplotlib.pyplot as plt
from uvlf_hod import HODModel

# Create a unified model with all parameters
model = HODModel(
    z=6.0,           # Redshift
    eps0=0.1,        # Star formation efficiency
    Mc=10**11.5,     # Characteristic halo mass [M_sun]
    a=0.6,           # Low-mass slope
    b=0.35,          # High-mass slope
    sigma_UV=0.35,   # UV magnitude scatter
    Mcut=10**10,     # Satellite cutoff mass [M_sun]
    Msat=10**12.5,   # Satellite normalization mass [M_sun]
    asat=1.0         # Satellite power-law slope
)

# Compute luminosity function
MUV = np.linspace(-22, -16, 20)
phi = model.luminosity_function(MUV)

# Plot
plt.semilogy(MUV, phi)
plt.xlabel('$M_{UV}$')
plt.ylabel('$\\Phi$ [Mpc$^{-3}$ mag$^{-1}$]')
plt.show()
```

### Galaxy Bias Calculation

```python
# Use the same model to compute bias
bias = model.galaxy_bias(MUV)

plt.plot(MUV, bias)
plt.xlabel('$M_{UV}$')
plt.ylabel('Galaxy Bias $b_g$')
plt.show()
```

### UV-Halo Mass Relation

```python
# UVHMR methods are inherited
Mh = 1e11  # M_sun
MUV = model.MUV(Mh)
sfr = model.sfr(Mh)

print(f"Halo mass {Mh:.2e} M_sun:")
print(f"  M_UV = {MUV:.2f}")
print(f"  SFR = {sfr:.1f} M_sun/yr")

# Inverse relation
Mh_recovered = model.Mhalo(MUV)
```

### Mean Properties

```python
# Compute mean properties for galaxies above a threshold
MUV_thresh = -20

mean_mass = model.mean_halo_mass(MUV_thresh)
mean_bias = model.mean_bias(MUV_thresh)

print(f"Galaxies brighter than {MUV_thresh}:")
print(f"  Mean halo mass: {10**mean_mass:.2e} M_sun")
print(f"  Mean bias: {mean_bias:.2f}")
```

### Updating Parameters

```python
# Dynamically update parameters
model.update_parameters(z=7.0, eps0=0.15)

# Recompute with new parameters
phi_new = model.luminosity_function(MUV)
```

## Documentation

Full documentation is available at [link to docs].

### Package Structure

```
uvlf_hod/
├── __init__.py          # Public API
├── config.py            # Configuration and constants
├── cosmology.py         # Halo mass function and bias
├── model.py             # Unified UVHMR and HOD models
├── luminosity.py        # UV luminosity and dust
└── models/
    └── parametrization.py  # Redshift parametrizations
```

### Key Classes

- **`HODModel`**: Main class combining UVHMR + HOD (recommended for most users)
- **`UVHMRModel`**: Base class for UV-halo mass relations only
- **`HaloMassFunction`**: Halo mass function with caching
- **`CosmologyConfig`**: Cosmology configuration

## Model Architecture

The package uses a clean inheritance hierarchy:

```
UVHMRModel (base class)
├── Handles UV-halo mass relations
├── Methods: sfr(), MUV(), Mhalo()
└── Parameters: z, eps0, Mc, a, b

HODModel (extends UVHMRModel)
├── Inherits all UVHMR methods
├── Adds occupation distributions
├── Methods: Ncen(), Nsat(), luminosity_function(), galaxy_bias()
└── Additional parameters: sigma_UV, Mcut, Msat, asat
```

## Examples

See the `examples/` directory for detailed examples:

- `basic_usage.ipynb`: Interactive Jupyter notebook with complete workflow
- `unified_model_example.py`: Comprehensive Python script
- `parameter_exploration.py`: Parameter sensitivity analysis
- `multi_redshift.py`: Redshift evolution studies

## Physics Background

### UV-Halo Mass Relation

The package connects halo mass to UV luminosity through:

1. **Star Formation**: SFR = ε(M_h) × f_b × dM_h/dt
2. **UV Luminosity**: L_UV = SFR / c_UV
3. **Dust Attenuation**: Self-consistent treatment using UV slope β

### HOD Model

Galaxy occupation follows:

- **Centrals**: N_cen(M_h | M_UV) via scatter around mean UVHMR
- **Satellites**: N_sat(M_h | M_UV) with cutoff mass and power-law

### Key Parameters

#### UVHMR Parameters

- **eps0**: Star formation efficiency normalization (typical: 0.05-0.2)
- **Mc**: Characteristic halo mass (typical: 10^11 - 10^12 M_sun)
- **a**: Low-mass slope (typical: 0.4-0.8)
- **b**: High-mass slope (typical: 0.2-0.5)

#### HOD Parameters

- **sigma_UV**: UV magnitude scatter (typical: 0.2-0.5 mag)
- **Mcut**: Cutoff mass for satellites (typical: 10^9 - 10^11 M_sun)
- **Msat**: Satellite normalization mass (typical: 10^11 - 10^13 M_sun)
- **asat**: Satellite power-law slope (typical: 0.8-1.2)

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
- Sabti et al. 2022
- Muñoz et al. 2023
- Bouwens et al. 2013, 2014
- And others...

## Contact

For questions or issues, please open an issue on GitHub or contact marko.shuntov@nbu.ku.dk.
