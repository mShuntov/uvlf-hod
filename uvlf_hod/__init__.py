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

# Convenience functions (backwards compatibility)
def compute_luminosity_function(MUV_arr, z, uvhmr_params, nsat_params, 
                                sigma_UV, add_dust=True):
    """Compute UV luminosity function.
    
    Convenience function that creates a HODModel and computes the LF.
    For better performance with multiple calculations, create a HODModel
    object directly.
    
    Parameters
    ----------
    MUV_arr : array_like
        Array of UV absolute magnitudes
    z : float
        Redshift
    uvhmr_params : list
        [eps0, log10(Mc), a, b]
    nsat_params : list
        [log10(Mcut), log10(Msat), asat]
    sigma_UV : float
        UV magnitude scatter
    add_dust : bool, optional
        Include dust attenuation
        
    Returns
    -------
    phi : ndarray
        Luminosity function Φ(M_UV) in Mpc^-3 mag^-1
    """
    model = HODModel(
        z=z,
        eps0=uvhmr_params[0],
        Mc=10**uvhmr_params[1],
        a=uvhmr_params[2],
        b=uvhmr_params[3],
        sigma_UV=sigma_UV,
        Mcut=10**nsat_params[0],
        Msat=10**nsat_params[1],
        asat=nsat_params[2],
        add_dust=add_dust
    )
    return model.luminosity_function(MUV_arr)


def compute_galaxy_bias(MUV_arr, z, uvhmr_params, nsat_params, 
                       sigma_UV, add_dust=True):
    """Compute galaxy bias.
    
    Convenience function that creates a HODModel and computes bias.
    
    Parameters
    ----------
    MUV_arr : array_like
        Array of UV absolute magnitudes
    z : float
        Redshift
    uvhmr_params : list
        [eps0, log10(Mc), a, b]
    nsat_params : list
        [log10(Mcut), log10(Msat), asat]
    sigma_UV : float
        UV magnitude scatter
    add_dust : bool, optional
        Include dust attenuation
        
    Returns
    -------
    bias : ndarray
        Galaxy bias
    """
    model = HODModel(
        z=z,
        eps0=uvhmr_params[0],
        Mc=10**uvhmr_params[1],
        a=uvhmr_params[2],
        b=uvhmr_params[3],
        sigma_UV=sigma_UV,
        Mcut=10**nsat_params[0],
        Msat=10**nsat_params[1],
        asat=nsat_params[2],
        add_dust=add_dust
    )
    return model.galaxy_bias(MUV_arr)


def compute_mean_halo_mass(MUV_thresh, z, uvhmr_params, nsat_params,
                           sigma_UV, add_dust=True):
    """Compute mean halo mass.
    
    Convenience function that creates a HODModel and computes mean mass.
    
    Parameters
    ----------
    MUV_thresh : float
        UV magnitude threshold
    z : float
        Redshift
    uvhmr_params : list
        [eps0, log10(Mc), a, b]
    nsat_params : list
        [log10(Mcut), log10(Msat), asat]
    sigma_UV : float
        UV magnitude scatter
    add_dust : bool, optional
        Include dust attenuation
        
    Returns
    -------
    log10_Mh : float
        Mean log10 halo mass
    """
    model = HODModel(
        z=z,
        eps0=uvhmr_params[0],
        Mc=10**uvhmr_params[1],
        a=uvhmr_params[2],
        b=uvhmr_params[3],
        sigma_UV=sigma_UV,
        Mcut=10**nsat_params[0],
        Msat=10**nsat_params[1],
        asat=nsat_params[2],
        add_dust=add_dust
    )
    return model.mean_halo_mass(MUV_thresh)


def compute_mean_bias(MUV_thresh, z, uvhmr_params, nsat_params,
                     sigma_UV, add_dust=True):
    """Compute mean galaxy bias.
    
    Convenience function that creates a HODModel and computes mean bias.
    
    Parameters
    ----------
    MUV_thresh : float
        UV magnitude threshold
    z : float
        Redshift
    uvhmr_params : list
        [eps0, log10(Mc), a, b]
    nsat_params : list
        [log10(Mcut), log10(Msat), asat]
    sigma_UV : float
        UV magnitude scatter
    add_dust : bool, optional
        Include dust attenuation
        
    Returns
    -------
    bias : float
        Mean galaxy bias
    """
    model = HODModel(
        z=z,
        eps0=uvhmr_params[0],
        Mc=10**uvhmr_params[1],
        a=uvhmr_params[2],
        b=uvhmr_params[3],
        sigma_UV=sigma_UV,
        Mcut=10**nsat_params[0],
        Msat=10**nsat_params[1],
        asat=nsat_params[2],
        add_dust=add_dust
    )
    return model.mean_bias(MUV_thresh)


__all__ = [
    # Version
    '__version__',
    
    # Main model classes
    'UVHMRModel',
    'HODModel',
    'UVHaloMassRelation',  # Backwards compatibility
    'UVLuminosityFunction',  # Backwards compatibility
    
    # Convenience functions
    'compute_luminosity_function',
    'compute_galaxy_bias',
    'compute_mean_halo_mass',
    'compute_mean_bias',
    
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