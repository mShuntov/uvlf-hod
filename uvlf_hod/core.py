"""Core physics calculations for UVHMR and HOD models."""

import numpy as np
from scipy import interpolate, integrate
from scipy.special import erf

from .config import BARYON_FRACTION
from .luminosity import MUV_from_SFR


def star_formation_efficiency(Mh, eps0, Mc, a, b):
    """Parametric form of instantaneous star formation efficiency.
    
    Double power-law form centered at characteristic mass Mc.
    
    Parameters
    ----------
    Mh : float or array_like
        Halo mass in M_sun
    eps0 : float
        Normalization (peak efficiency)
    Mc : float
        Characteristic halo mass in M_sun
    a : float
        Low-mass power-law slope
    b : float
        High-mass power-law slope
        
    Returns
    -------
    epsilon : float or ndarray
        Star formation efficiency
    """
    Mh = np.atleast_1d(Mh)
    return 2 * eps0 / ((Mh / Mc)**(-a) + (Mh / Mc)**b)


def halo_accretion_rate(Mh, z):
    """Halo mass accretion rate.
    
    Parameters
    ----------
    Mh : float or array_like
        Halo mass in M_sun
    z : float
        Redshift
        
    Returns
    -------
    dMh_dt : float or ndarray
        Halo mass accretion rate in M_sun/yr
    """
    Mh = np.atleast_1d(Mh)
    return 0.03e-9 * Mh * (Mh / 1e12)**0.14 * (1 + z)**2.5


def sfr_from_halo_mass(Mh, z, eps0, Mc, a, b):
    """Star formation rate from halo mass.
    
    SFR = epsilon * f_b * dMh/dt
    
    Parameters
    ----------
    Mh : float or array_like
        Halo mass in M_sun
    z : float
        Redshift
    eps0, Mc, a, b : float
        Star formation efficiency parameters
        
    Returns
    -------
    sfr : float or ndarray
        Star formation rate in M_sun/yr
    """
    f_star = star_formation_efficiency(Mh, eps0, Mc, a, b)
    return BARYON_FRACTION * f_star * halo_accretion_rate(Mh, z)


def MUV_from_halo_mass(Mh, z, eps0, Mc, a, b, add_dust=True):
    """UV absolute magnitude from halo mass (UVHMR).
    
    Parameters
    ----------
    Mh : float or array_like
        Halo mass in M_sun
    z : float
        Redshift
    eps0, Mc, a, b : float
        Star formation efficiency parameters
    add_dust : bool, optional
        Whether to add dust attenuation
        
    Returns
    -------
    MUV : float or ndarray
        UV absolute magnitude
    """
    sfr = sfr_from_halo_mass(Mh, z, eps0, Mc, a, b)
    return MUV_from_SFR(sfr, z, add_dust=add_dust)


def halo_mass_from_MUV(MUV, z, eps0, Mc, a, b, add_dust=True,
                       M_min=1e6, M_max=1e16, num_points=1000):
    """Inverse UVHMR: halo mass from UV magnitude.
    
    Parameters
    ----------
    MUV : float or array_like
        UV absolute magnitude
    z : float
        Redshift
    eps0, Mc, a, b : float
        Star formation efficiency parameters
    add_dust : bool, optional
        Whether dust was included in MUV
    M_min, M_max : float, optional
        Halo mass range for interpolation
    num_points : int, optional
        Number of points for interpolation
        
    Returns
    -------
    Mh : float or ndarray
        Halo mass in M_sun
    """
    # Create forward relation
    Mhalo = np.logspace(np.log10(M_min), np.log10(M_max), num_points)
    MUV_forward = MUV_from_halo_mass(Mhalo, z, eps0, Mc, a, b, add_dust)
    
    # Invert using interpolation
    interp_func = interpolate.interp1d(MUV_forward, Mhalo, 
                                       fill_value='extrapolate')
    
    return interp_func(MUV)


def P_MUV_given_Mhalo(MUV, MUV_mean, sigma_UV):
    """Probability distribution of MUV given halo mass.
    
    Gaussian scatter around mean MUV(Mhalo) relation.
    
    Parameters
    ----------
    MUV : float or array_like
        UV absolute magnitude
    MUV_mean : float
        Mean MUV for given halo mass
    sigma_UV : float
        Scatter in magnitudes
        
    Returns
    -------
    P : float or ndarray
        Probability density
    """
    return (1.0 / (np.sqrt(2 * np.pi) * sigma_UV) *
            np.exp(-0.5 * ((MUV - MUV_mean) / sigma_UV)**2))


class UVHaloMassRelation:
    """Class for handling UV-halo mass relation calculations.
    
    Attributes
    ----------
    z : float
        Redshift
    eps0, Mc, a, b : float
        UVHMR parameters
    add_dust : bool
        Whether to include dust
    """
    
    def __init__(self, z, eps0, Mc, a, b, add_dust=True):
        """Initialize UVHMR.
        
        Parameters
        ----------
        z : float
            Redshift
        eps0 : float
            Star formation efficiency normalization
        Mc : float
            Characteristic halo mass
        a : float
            Low-mass slope
        b : float
            High-mass slope
        add_dust : bool, optional
            Include dust attenuation
        """
        self.z = z
        self.eps0 = eps0
        self.Mc = Mc
        self.a = a
        self.b = b
        self.add_dust = add_dust
    
    def MUV(self, Mh):
        """Get MUV for halo mass."""
        return MUV_from_halo_mass(Mh, self.z, self.eps0, self.Mc, 
                                   self.a, self.b, self.add_dust)
    
    def Mhalo(self, MUV):
        """Get halo mass for MUV."""
        return halo_mass_from_MUV(MUV, self.z, self.eps0, self.Mc,
                                   self.a, self.b, self.add_dust)
    
    def sfr(self, Mh):
        """Get SFR for halo mass."""
        return sfr_from_halo_mass(Mh, self.z, self.eps0, self.Mc,
                                   self.a, self.b)
    
    def __repr__(self):
        return (f"UVHaloMassRelation(z={self.z}, eps0={self.eps0}, "
                f"Mc={self.Mc:.2e}, a={self.a}, b={self.b})")
