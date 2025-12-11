"""Cosmology-related functions for halo mass function and bias calculations."""

import numpy as np
from colossus.lss import mass_function, bias

from .config import HMF_NUM_POINTS


def get_halo_mass_function(z, M_min=6, M_max=15, num_points=None, 
                           mdef='vir', model='watson13'):
    """Compute the halo mass function at redshift z.
    
    Parameters
    ----------
    z : float
        Redshift
    M_min : float, optional
        Minimum log10(M/M_sun). Default is 6.
    M_max : float, optional
        Maximum log10(M/M_sun). Default is 15.
    num_points : int, optional
        Number of mass bins. Default is from config.
    mdef : str, optional
        Mass definition: 'vir', 'fof', '200m', etc. Default is 'vir'.
    model : str, optional
        HMF model: 'watson13', 'tinker08', 'reed07', etc. Default is 'watson13'.
        
    Returns
    -------
    log10_mass : ndarray
        Array of log10(M/M_sun) including h factor
    hmf : ndarray
        Halo mass function dn/dln(M) in units of h^3 Mpc^-3
        
    Notes
    -----
    The mass is returned in units that include h: M [M_sun/h]
    The HMF is returned in units that include h^3: dn/dlnM [h^3 Mpc^-3]
    """
    if num_points is None:
        num_points = HMF_NUM_POINTS
    
    # Create mass array (without h initially)
    log10_mass_noh = np.linspace(M_min, M_max, num_points)
    mass_noh = 10**log10_mass_noh
    
    # Get HMF from Colossus (in comoving units without h)
    hmf_noh = mass_function.massFunction(
        mass_noh, z, 
        mdef=mdef, 
        model=model, 
        q_out='dndlnM'
    )
    
    # Convert to include h factors
    h = 0.7  # H0/100
    log10_mass = log10_mass_noh - np.log10(h)  # M [M_sun/h]
    hmf = hmf_noh * np.log(10) * h**3  # dn/dlnM [h^3 Mpc^-3]
    
    return log10_mass, hmf


def get_halo_bias(M, z, mdef='vir', model='tinker10'):
    """Compute halo bias as a function of mass and redshift.
    
    Parameters
    ----------
    M : float or array_like
        Halo mass in M_sun (without h factor)
    z : float
        Redshift
    mdef : str, optional
        Mass definition. Default is 'vir'.
    model : str, optional
        Bias model. Default is 'tinker10'.
        
    Returns
    -------
    b : float or ndarray
        Halo bias
    """
    return bias.haloBias(M, model=model, z=z, mdef=mdef)


class HaloMassFunction:
    """Class for computing and caching halo mass functions.
    
    This class provides a convenient interface for HMF calculations
    with optional caching of results.
    
    Attributes
    ----------
    z : float
        Redshift
    M_min : float
        Minimum log10(M/M_sun)
    M_max : float
        Maximum log10(M/M_sun)
    num_points : int
        Number of mass bins
    mdef : str
        Mass definition
    model : str
        HMF model
    """
    
    def __init__(self, z, M_min=6, M_max=15, num_points=None,
                 mdef='vir', model='watson13'):
        """Initialize HMF calculator.
        
        Parameters
        ----------
        z : float
            Redshift
        M_min : float, optional
            Minimum log10(M/M_sun)
        M_max : float, optional
            Maximum log10(M/M_sun)
        num_points : int, optional
            Number of mass bins
        mdef : str, optional
            Mass definition
        model : str, optional
            HMF model
        """
        self.z = z
        self.M_min = M_min
        self.M_max = M_max
        self.num_points = num_points if num_points else HMF_NUM_POINTS
        self.mdef = mdef
        self.model = model
        
        # Cache
        self._log10_mass = None
        self._hmf = None
        self._bias = None
    
    @property
    def log10_mass(self):
        """Get log10 mass array."""
        if self._log10_mass is None:
            self._compute()
        return self._log10_mass
    
    @property
    def mass(self):
        """Get linear mass array."""
        return 10**self.log10_mass
    
    @property
    def hmf(self):
        """Get halo mass function."""
        if self._hmf is None:
            self._compute()
        return self._hmf
    
    @property
    def halo_bias(self):
        """Get halo bias."""
        if self._bias is None:
            self._bias = get_halo_bias(self.mass, self.z, 
                                       mdef=self.mdef, model='tinker10')
        return self._bias
    
    def _compute(self):
        """Compute HMF and cache results."""
        self._log10_mass, self._hmf = get_halo_mass_function(
            self.z, self.M_min, self.M_max, self.num_points,
            self.mdef, self.model
        )
    
    def __repr__(self):
        return (f"HaloMassFunction(z={self.z}, M_min={self.M_min}, "
                f"M_max={self.M_max}, model='{self.model}')")
