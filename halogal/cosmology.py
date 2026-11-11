"""Cosmology-related functions for halo mass function and bias calculations."""

import numpy as np
from colossus.cosmology import cosmology as colossus_cosmology
from colossus.lss import mass_function, bias

from .config import HMF_NUM_POINTS

# Colossus returns dn/dlnM; every integral in this package uses log10(M) as the
# integration variable, so the package-wide convention for an "hmf" array is
# dn/dlog10M [Mpc^-3 dex^-1]. The two differ by exactly this factor.
LN10 = np.log(10)


def current_h():
    """Return h = H0/100 of the currently active Colossus cosmology.

    Returns
    -------
    h : float
        Dimensionless Hubble parameter.

    Notes
    -----
    The M_sun <-> M_sun/h and Mpc^-3 <-> (h/Mpc)^3 conversions below must use
    the same h that Colossus evaluated the HMF with. Hardcoding h=0.7 left the
    mass grid bit-identical when the user changed H0 via ``CosmologyConfig``,
    which silently biased masses and amplitudes.
    """
    return colossus_cosmology.getCurrent().h


def get_halo_mass_function(z, M_min=6, M_max=15, num_points=None,
                           mdef='vir', model='watson13', h=None):
    """Compute the halo mass function at redshift z.

    Parameters
    ----------
    z : float
        Redshift
    M_min : float, optional
        Minimum log10(M/M_sun) in physical units. Default is 6.
    M_max : float, optional
        Maximum log10(M/M_sun) in physical units. Default is 15.
    num_points : int, optional
        Number of mass bins. Default is from config.
    mdef : str, optional
        Mass definition: 'vir', 'fof', '200m', etc. Default is 'vir'.
    model : str, optional
        HMF model: 'watson13', 'tinker08', 'reed07', etc. Default is 'watson13'.
    h : float, optional
        Dimensionless Hubble parameter used for the unit conversions. Default
        is taken from the active Colossus cosmology (see `current_h`).

    Returns
    -------
    log10_mass : ndarray
        Array of log10(M/M_sun) in physical units
    hmf : ndarray
        Halo mass function dn/dlog10(M) in physical units of Mpc^-3 dex^-1

    Notes
    -----
    The returned quantity is dn/dlog10M, NOT dn/dlnM: it is meant to be
    integrated directly against ``log10_mass``. Earlier versions of this
    docstring said dn/dlnM while the code already returned dn/dlog10M, which is
    a factor ln(10) = 2.3026 (0.362 dex) trap for external callers.

    Colossus expects mass in M_sun/h and returns dn/dlnM in (h/Mpc)^3:
    - Mass: M [M_sun/h] = M [M_sun] * h  (factor out h by multiplying)
    - HMF: dn/dlnM [Mpc^-3] = dn/dlnM [(h/Mpc)^3] * h^3  (evaluate h)
    - Measure: dn/dlog10M = dn/dlnM * ln(10)
    """
    if num_points is None:
        num_points = HMF_NUM_POINTS

    if h is None:
        h = current_h()

    # Create mass array in physical units [M_sun]
    log10_mass_physical = np.linspace(M_min, M_max, num_points)
    mass_physical = 10**log10_mass_physical
    
    # Convert to comoving units for Colossus [M_sun/h]
    # To factor out h: M [M_sun/h] = M [M_sun] * h
    mass_comoving = mass_physical * h
    
    # Get HMF from Colossus in comoving units [(h/Mpc)^3]
    hmf_comoving = mass_function.massFunction(
        mass_comoving, z, 
        mdef=mdef, 
        model=model, 
        q_out='dndlnM'
    )
    
    # Evaluate h to get Mpc^-3, then change measure dlnM -> dlog10M so the
    # result can be integrated against log10_mass_physical directly
    hmf_physical = hmf_comoving * LN10 * h**3

    return log10_mass_physical, hmf_physical


def get_halo_bias(M, z, mdef='vir', model='tinker10', h=None):
    """Compute halo bias as a function of mass and redshift.
    
    Parameters
    ----------
    M : float or array_like
        Halo mass in M_sun (physical units)
    z : float
        Redshift
    mdef : str, optional
        Mass definition. Default is 'vir'.
    model : str, optional
        Bias model. Default is 'tinker10'.
    h : float, optional
        Dimensionless Hubble parameter used for the unit conversion. Default
        is taken from the active Colossus cosmology (see `current_h`).

    Returns
    -------
    b : float or ndarray
        Halo bias (dimensionless)

    Notes
    -----
    Colossus expects mass in M_sun/h.
    To factor out h: M [M_sun/h] = M [M_sun] * h
    """
    if h is None:
        h = current_h()
    M_comoving = np.atleast_1d(M) * h  # Convert to M_sun/h
    return bias.haloBias(M_comoving, model=model, z=z, mdef=mdef)


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