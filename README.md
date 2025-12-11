# uvlf-hod

[![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub](https://img.shields.io/github/stars/mshuntov/uvlf-hod?style=social)](https://github.com/mshuntov/uvlf-hod)

A Python package for modeling UV luminosity functions and galaxy clustering using Halo Occupation Distribution (HOD) models for high-redshift galaxies.

## Features

- **UV-Halo Mass Relation (UVHMR)**: Connect halo mass to UV luminosity through star formation
- **Halo Occupation Distribution**: Model central and satellite galaxy populations
- **Luminosity Functions**: Compute UV luminosity functions at high redshift
- **Galaxy Bias**: Calculate galaxy clustering bias
- **Dust Attenuation**: Self-consistent treatment following Bouwens+2013-14
- **Flexible Parametrization**: Redshift-dependent model parameters
- **Modern Python**: Clean API with type hints and comprehensive documentation

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

### Basic UV Luminosity Function

```python
import numpy as np
from uvlf_hod import UVLuminosityFunction

# Define redshift and magnitude range
z = 6.0
MUV = np.linspace(-22, -16, 20)

# Set up model parameters
uvhmr_params = [0.1, 11.5, 0.6, 0.35]  # [eps0, log10(Mc), a, b]
nsat_params = [10.0, 12.5, 1.0]  # [log10(Mcut), log10(Msat), asat]
sigma_UV = 0.35

# Calculate luminosity function
uvlf = UVLuminosityFunction(z, uvhmr_params, nsat_params, sigma_UV)
phi = uvlf.compute(MUV)

# Plot
import matplotlib.pyplot as plt
plt.semilogy(MUV, phi)
plt.xlabel('$M_{UV}$')
plt.ylabel('$\\Phi$ [Mpc$^{-3}$ mag$^{-1}$]')
plt.show()
```

### Galaxy Bias Calculation

```python
from uvlf_hod import compute_galaxy_bias

# Calculate bias for magnitude bins
bias = compute_galaxy_bias(MUV, z, uvhmr_params, nsat_params, sigma_UV)

plt.plot(MUV, bias)
plt.xlabel('$M_{UV}$')
plt.ylabel('Galaxy Bias $b_g$')
plt.show()
```

### UV-Halo Mass Relation

```python
from uvlf_hod.core import UVHaloMassRelation

# Create UVHMR instance
uvhmr = UVHaloMassRelation(z=6.0, eps0=0.1, Mc=10**11.5, a=0.6, b=0.35)

# Get UV magnitude for halo mass
Mh = 1e11  # M_sun
MUV = uvhmr.MUV(Mh)
print(f"Halo mass {Mh:.2e} M_sun -> MUV = {MUV:.2f}")

# Inverse relation
Mh_inv = uvhmr.Mhalo(MUV)
print(f"MUV = {MUV:.2f} -> Halo mass {Mh_inv:.2e} M_sun")
```

## Documentation

Full documentation is available at [link to docs].

### Package Structure

```
uvlf_hod/
├── __init__.py          # Public API
├── config.py            # Configuration and constants
├── cosmology.py         # Halo mass function and bias
├── core.py              # Core UVHMR calculations
├── hod.py               # HOD models and LF calculations
├── luminosity.py        # UV luminosity and dust
└── models/
    ├── parametrization.py  # Redshift parametrizations
    └── occupation.py       # Occupation function models
```

### Key Classes

- `UVLuminosityFunction`: Main class for computing luminosity functions
- `UVHaloMassRelation`: UV-halo mass relation calculations
- `HaloMassFunction`: Halo mass function with caching
- `HODModel`: Base class for HOD models
- `CosmologyConfig`: Cosmology configuration

## Examples

See the `examples/` directory for detailed examples:

- `basic_usage.py`: Basic luminosity function calculation
- `compute_luminosity_function.py`: Full LF with error handling
- `galaxy_bias_calculation.py`: Computing galaxy bias
- `redshift_evolution.py`: Parameter evolution with redshift
- `fitting_observations.py`: Fitting to observational data

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

- `eps0`: Star formation efficiency normalization
- `Mc`: Characteristic halo mass
- `a, b`: UVHMR slopes (low/high mass)
- `sigma_UV`: UV magnitude scatter
- `Mcut, Msat, asat`: Satellite parameters

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
