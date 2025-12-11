"""
uvlf-hod: UV Luminosity Functions with Halo Occupation Distribution models.

This package provides tools for modeling high-redshift galaxy populations
using UV-halo mass relations and HOD models.

Main components:
    - UV-Halo Mass Relation (UVHMR)
    - Halo Occupation Distribution (HOD) models
    - UV Luminosity Function calculations
    - Galaxy bias calculations
    - Dust attenuation models
"""

from .__version__ import __version__

# Core functionality
from .core import (
    UVHaloMassRelation,
    star_formation_efficiency,
    halo_accretion_rate,
    sfr_from_halo_mass,
    MUV_from_halo_mass,
    halo_mass_from_MUV,
)

# Cosmology
from .cosmology import (
    HaloMassFunction,
    get_halo_mass_function,
    get_halo_bias,
)

# Luminosity
from .luminosity import (
    beta_color,
    dust_attenuation,
    MUV_from_SFR,
    SFR_from_MUV,
)

# HOD and LF calculations
from .hod import (
    UVLuminosityFunction,
    compute_luminosity_function,
    compute_galaxy_bias,
    compute_mean_halo_mass,
    compute_mean_bias,
    occupation_central,
    occupation_satellite,
)

# Configuration
from .config import (
    CosmologyConfig,
    default_cosmology,
    DEFAULT_HOD_PARAMS,
)

# Models
from .models.parametrization import (
    RedshiftParametrization,
    eps0_fz,
    Mc_fz,
    a_fz,
    b_fz,
    sigma_UV_fz,
    Mcut_fz,
    Msat_fz,
    asat_fz,
)


__all__ = [
    # Version
    '__version__',
    
    # Core
    'UVHaloMassRelation',
    'star_formation_efficiency',
    'halo_accretion_rate',
    'sfr_from_halo_mass',
    'MUV_from_halo_mass',
    'halo_mass_from_MUV',
    
    # Cosmology
    'HaloMassFunction',
    'get_halo_mass_function',
    'get_halo_bias',
    
    # Luminosity
    'beta_color',
    'dust_attenuation',
    'MUV_from_SFR',
    'SFR_from_MUV',
    
    # HOD
    'UVLuminosityFunction',
    'compute_luminosity_function',
    'compute_galaxy_bias',
    'compute_mean_halo_mass',
    'compute_mean_bias',
    'occupation_central',
    'occupation_satellite',
    
    # Config
    'CosmologyConfig',
    'default_cosmology',
    'DEFAULT_HOD_PARAMS',
    
    # Models
    'RedshiftParametrization',
    'eps0_fz',
    'Mc_fz',
    'a_fz',
    'b_fz',
    'sigma_UV_fz',
    'Mcut_fz',
    'Msat_fz',
    'asat_fz',
]


# Package metadata
__author__ = "Your Name"
__email__ = "your.email@example.com"
__url__ = "https://github.com/yourusername/uvlf-hod"
__description__ = "UV Luminosity Functions with HOD models for high-redshift galaxies"