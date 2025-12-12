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

# Main model classes (new unified interface)
from .model import (
    UVHMRModel,
    HODModel,
    UVHaloMassRelation,  # Alias for backwards compatibility
    UVLuminosityFunction,  # Alias for backwards compatibility
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

# Configuration
from .config import (
    CosmologyConfig,
    default_cosmology,
    DEFAULT_HOD_PARAMS,
    DEFAULT_REDSHIFT_EVOLUTION,
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
    
    # Main model classes
    'UVHMRModel',
    'HODModel',
    
    # Cosmology
    'HaloMassFunction',
    'get_halo_mass_function',
    'get_halo_bias',
    
    # Luminosity
    'beta_color',
    'dust_attenuation',
    'MUV_from_SFR',
    'SFR_from_MUV',
    
    # Config
    'CosmologyConfig',
    'default_cosmology',
    'DEFAULT_HOD_PARAMS',
    'DEFAULT_REDSHIFT_EVOLUTION',
    
    # Redshift parametrization
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
__author__ = "Marko Shuntov"
__email__ = "marko.shuntov@nbi.ku.dk"
__url__ = "https://github.com/mshuntov/uvlf-hod"
__description__ = "UV Luminosity Functions with HOD models for high-redshift galaxies"