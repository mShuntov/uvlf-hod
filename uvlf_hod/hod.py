"""Halo Occupation Distribution models and luminosity function calculations."""

import numpy as np
from scipy import integrate
from scipy.special import erf

from .core import MUV_from_halo_mass
from .cosmology import get_halo_mass_function, get_halo_bias


def occupation_central(Mh, MUV_thresh, z, eps0, Mc, a, b, sigma_UV, add_dust=True):
    """Central galaxy occupation function.
    
    Average number of central galaxies with M_UV < MUV_thresh in halo of mass Mh.
    
    Parameters
    ----------
    Mh : float or array_like
        Halo mass in M_sun
    MUV_thresh : float
        UV magnitude threshold (brighter than)
    z : float
        Redshift
    eps0, Mc, a, b : float
        UVHMR parameters
    sigma_UV : float
        UV magnitude scatter
    add_dust : bool, optional
        Include dust attenuation
        
    Returns
    -------
    N_cen : float or ndarray
        Average number of central galaxies
    """
    MUV_mean = MUV_from_halo_mass(Mh, z, eps0, Mc, a, b, add_dust)
    return 0.5 * (1 - erf((MUV_mean - MUV_thresh) / (np.sqrt(2) * sigma_UV)))


def occupation_satellite(Mh, MUV_thresh, z, eps0, Mc, a, b, sigma_UV,
                         Mcut, Msat, asat, add_dust=True):
    """Satellite galaxy occupation function.
    
    Average number of satellite galaxies with M_UV < MUV_thresh in halo of mass Mh.
    
    Parameters
    ----------
    Mh : float or array_like
        Halo mass in M_sun
    MUV_thresh : float
        UV magnitude threshold
    z : float
        Redshift
    eps0, Mc, a, b : float
        UVHMR parameters
    sigma_UV : float
        UV magnitude scatter
    Mcut : float
        Cutoff mass for satellites
    Msat : float
        Satellite normalization mass
    asat : float
        Satellite power-law slope
    add_dust : bool, optional
        Include dust attenuation
        
    Returns
    -------
    N_sat : float or ndarray
        Average number of satellite galaxies
    """
    Mh = np.atleast_1d(Mh)
    N_cen = occupation_central(Mh, MUV_thresh, z, eps0, Mc, a, b, sigma_UV, add_dust)
    
    # Satellites only exist above cutoff mass
    N_sat = np.zeros_like(Mh)
    mask = Mh > Mcut
    N_sat[mask] = ((Mh[mask] - Mcut) / Msat)**asat
    
    return N_sat * N_cen


def _galaxy_number_density_integrand(log10_Mh, MUV, dMUV, z, uvhmr_params,
                                     nsat_params, sigma_UV, add_dust, hmf):
    """Calculate integrand for galaxy number density."""
    eps0, Mc, a, b = uvhmr_params
    Mcut, Msat, asat = nsat_params
    Mh = 10**log10_Mh
    
    # Number in magnitude bin [MUV - dMUV/2, MUV + dMUV/2]
    N_cen_bin = (occupation_central(Mh, MUV + dMUV/2, z, eps0, Mc, a, b, sigma_UV, add_dust) -
                 occupation_central(Mh, MUV - dMUV/2, z, eps0, Mc, a, b, sigma_UV, add_dust))
    
    N_sat_bin = (occupation_satellite(Mh, MUV + dMUV/2, z, eps0, Mc, a, b, sigma_UV,
                                      10**Mcut, 10**Msat, asat, add_dust) -
                 occupation_satellite(Mh, MUV - dMUV/2, z, eps0, Mc, a, b, sigma_UV,
                                      10**Mcut, 10**Msat, asat, add_dust))
    
    return (N_cen_bin + N_sat_bin) * hmf


def compute_number_density_at_MUV(MUV, z, uvhmr_params, nsat_params, sigma_UV,
                                   add_dust=True, dMUV=0.05):
    """Compute galaxy number density at specific MUV.
    
    Parameters
    ----------
    MUV : float
        UV absolute magnitude
    z : float
        Redshift
    uvhmr_params : list
        [eps0, log10(Mc), a, b]
    nsat_params : list
        [log10(Mcut), log10(Msat), asat]
    sigma_UV : float
        UV magnitude scatter
    add_dust : bool, optional
        Include dust
    dMUV : float, optional
        Magnitude bin width
        
    Returns
    -------
    dn_dMUV : float
        Number density per magnitude
    """
    # Get halo mass function
    log10_Mh, hmf = get_halo_mass_function(z, M_min=6, M_max=15, num_points=2048)
    
    # Unpack parameters
    eps0, log_Mc, a, b = uvhmr_params
    Mc = 10**log_Mc
    uvhmr_params_linear = [eps0, Mc, a, b]
    
    # Calculate integrand
    integrand = _galaxy_number_density_integrand(
        log10_Mh, MUV, dMUV, z, uvhmr_params_linear,
        nsat_params, sigma_UV, add_dust, hmf
    )
    
    # Integrate over halo mass
    result = integrate.simpson(integrand, x=log10_Mh) / dMUV
    
    return result


def compute_luminosity_function(MUV_arr, z, uvhmr_params, nsat_params, sigma_UV,
                                add_dust=True):
    """Compute UV luminosity function.
    
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
    Phi : ndarray
        Luminosity function Φ(M_UV) in Mpc^-3 mag^-1
        
    Examples
    --------
    >>> import numpy as np
    >>> MUV = np.linspace(-22, -16, 20)
    >>> uvhmr_params = [0.1, 11.5, 0.6, 0.35]
    >>> nsat_params = [10.0, 12.5, 1.0]
    >>> phi = compute_luminosity_function(MUV, 6.0, uvhmr_params, nsat_params, 0.35)
    """
    MUV_arr = np.atleast_1d(MUV_arr)
    
    phi = np.array([
        compute_number_density_at_MUV(muv, z, uvhmr_params, nsat_params, 
                                      sigma_UV, add_dust)
        for muv in MUV_arr
    ])
    
    return phi


def compute_galaxy_bias(MUV_arr, z, uvhmr_params, nsat_params, sigma_UV,
                       add_dust=True):
    """Compute galaxy bias as a function of UV magnitude.
    
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
        Galaxy bias as function of M_UV
        
    Examples
    --------
    >>> import numpy as np
    >>> MUV = np.linspace(-22, -16, 20)
    >>> uvhmr_params = [0.1, 11.5, 0.6, 0.35]
    >>> nsat_params = [10.0, 12.5, 1.0]
    >>> bias = compute_galaxy_bias(MUV, 6.0, uvhmr_params, nsat_params, 0.35)
    """
    # Unpack and prepare parameters
    eps0, log_Mc, a, b = uvhmr_params
    Mc = 10**log_Mc
    Mcut_log, Msat_log, asat = nsat_params
    
    # Get HMF and bias
    log10_Mh, hmf = get_halo_mass_function(z, M_min=6, M_max=15, num_points=2048)
    bias_halo = get_halo_bias(10**log10_Mh, z, mdef='vir', model='tinker10')
    
    dMUV = 0.05
    MUV_arr = np.atleast_1d(MUV_arr)
    bg = np.zeros_like(MUV_arr)
    
    for i, muv in enumerate(MUV_arr):
        # Occupation in bin
        N_cen_bin = (occupation_central(10**log10_Mh, muv + dMUV/2, z, 
                                        eps0, Mc, a, b, sigma_UV, add_dust) -
                     occupation_central(10**log10_Mh, muv - dMUV/2, z,
                                        eps0, Mc, a, b, sigma_UV, add_dust))
        
        N_sat_bin = (occupation_satellite(10**log10_Mh, muv + dMUV/2, z,
                                         eps0, Mc, a, b, sigma_UV,
                                         10**Mcut_log, 10**Msat_log, asat, add_dust) -
                     occupation_satellite(10**log10_Mh, muv - dMUV/2, z,
                                         eps0, Mc, a, b, sigma_UV,
                                         10**Mcut_log, 10**Msat_log, asat, add_dust))
        
        # Weighted bias
        bg_integrand = (N_cen_bin + N_sat_bin) * hmf * bias_halo
        ngal_integrand = (N_cen_bin + N_sat_bin) * hmf
        
        ngal = integrate.simpson(ngal_integrand, x=log10_Mh)
        
        if ngal > 0:
            bg[i] = integrate.simpson(bg_integrand, x=log10_Mh) / ngal
        else:
            bg[i] = 0.0
    
    return bg


def compute_mean_halo_mass(MUV_thresh, z, uvhmr_params, nsat_params, sigma_UV,
                           add_dust=True):
    """Compute mean halo mass for galaxies brighter than MUV threshold.
    
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
    log10_Mh_mean : float
        Mean log10 halo mass
    """
    # Unpack parameters
    eps0, log_Mc, a, b = uvhmr_params
    Mc = 10**log_Mc
    Mcut_log, Msat_log, asat = nsat_params
    
    # Get HMF
    log10_Mh, hmf = get_halo_mass_function(z, M_min=6, M_max=15, num_points=2048)
    
    # Occupation above threshold
    N_cen = occupation_central(10**log10_Mh, MUV_thresh, z, eps0, Mc, a, b, 
                              sigma_UV, add_dust)
    N_sat = occupation_satellite(10**log10_Mh, MUV_thresh, z, eps0, Mc, a, b,
                                sigma_UV, 10**Mcut_log, 10**Msat_log, asat, add_dust)
    
    # Mean mass
    mh_integrand = (N_cen + N_sat) * hmf * 10**log10_Mh
    ngal_integrand = (N_cen + N_sat) * hmf
    
    ngal = integrate.simpson(ngal_integrand, x=log10_Mh)
    
    if ngal > 0:
        mh_mean = integrate.simpson(mh_integrand, x=log10_Mh) / ngal
        return np.log10(mh_mean)
    else:
        return np.nan


def compute_mean_bias(MUV_thresh, z, uvhmr_params, nsat_params, sigma_UV,
                     add_dust=True):
    """Compute mean galaxy bias for galaxies brighter than MUV threshold.
    
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
    bias_mean : float
        Mean galaxy bias
    """
    # Unpack parameters
    eps0, log_Mc, a, b = uvhmr_params
    Mc = 10**log_Mc
    Mcut_log, Msat_log, asat = nsat_params
    
    # Get HMF and bias
    log10_Mh, hmf = get_halo_mass_function(z, M_min=6, M_max=15, num_points=2048)
    bias_halo = get_halo_bias(10**log10_Mh, z, mdef='vir', model='tinker10')
    
    # Occupation above threshold
    N_cen = occupation_central(10**log10_Mh, MUV_thresh, z, eps0, Mc, a, b,
                              sigma_UV, add_dust)
    N_sat = occupation_satellite(10**log10_Mh, MUV_thresh, z, eps0, Mc, a, b,
                                sigma_UV, 10**Mcut_log, 10**Msat_log, asat, add_dust)
    
    # Mean bias
    bg_integrand = (N_cen + N_sat) * hmf * bias_halo
    ngal_integrand = (N_cen + N_sat) * hmf
    
    ngal = integrate.simpson(ngal_integrand, x=log10_Mh)
    
    if ngal > 0:
        return integrate.simpson(bg_integrand, x=log10_Mh) / ngal
    else:
        return np.nan



class UVLuminosityFunction:
    """Main class for computing UV luminosity functions.
    
    This class provides a convenient interface for calculating UV luminosity
    functions with HOD models, including support for parameter variations and
    multi-redshift calculations.
    
    Parameters
    ----------
    z : float
        Redshift
    uvhmr_params : list
        UVHMR parameters [eps0, log10(Mc), a, b]
    nsat_params : list
        Satellite parameters [log10(Mcut), log10(Msat), asat]
    sigma_UV : float
        UV magnitude scatter
    add_dust : bool, optional
        Whether to include dust attenuation. Default is True.
        
    Attributes
    ----------
    z : float
        Redshift
    eps0, Mc, a, b : float
        UVHMR parameters (Mc is converted from log)
    Mcut, Msat, asat : float
        Satellite parameters (Mcut and Msat converted from log)
    sigma_UV : float
        UV magnitude scatter
    add_dust : bool
        Whether dust is included
        
    Examples
    --------
    >>> import numpy as np
    >>> from uvlf_hod import UVLuminosityFunction
    >>> 
    >>> # Set up model
    >>> uvhmr_params = [0.1, 11.5, 0.6, 0.35]
    >>> nsat_params = [10.0, 12.5, 1.0]
    >>> uvlf = UVLuminosityFunction(6.0, uvhmr_params, nsat_params, 0.35)
    >>> 
    >>> # Compute luminosity function
    >>> MUV = np.linspace(-22, -16, 20)
    >>> phi = uvlf.compute(MUV)
    >>> 
    >>> # Get galaxy bias
    >>> bias = uvlf.compute_bias(MUV)
    """
    
    def __init__(self, z, uvhmr_params, nsat_params, sigma_UV, add_dust=True):
        """Initialize UV luminosity function calculator."""
        self.z = z
        
        # Unpack and store UVHMR parameters
        self.eps0 = uvhmr_params[0]
        self.Mc = 10**uvhmr_params[1]  # Convert from log
        self.a = uvhmr_params[2]
        self.b = uvhmr_params[3]
        
        # Unpack and store satellite parameters
        self.Mcut = 10**nsat_params[0]  # Convert from log
        self.Msat = 10**nsat_params[1]  # Convert from log
        self.asat = nsat_params[2]
        
        # Store other parameters
        self.sigma_UV = sigma_UV
        self.add_dust = add_dust
        
        # Cache
        self._hmf = None
        self._log10_mass = None
    
    @property
    def uvhmr_params(self):
        """Get UVHMR parameters in standard format."""
        return [self.eps0, self.Mc, self.a, self.b]
    
    @property
    def nsat_params(self):
        """Get satellite parameters in log format."""
        return [np.log10(self.Mcut), np.log10(self.Msat), self.asat]
    
    def _get_hmf(self):
        """Get or compute halo mass function."""
        if self._hmf is None or self._log10_mass is None:
            self._log10_mass, self._hmf = get_halo_mass_function(
                self.z, M_min=6, M_max=15, num_points=2048
            )
        return self._log10_mass, self._hmf
    
    def compute(self, MUV_arr):
        """Compute UV luminosity function.
        
        Parameters
        ----------
        MUV_arr : array_like
            Array of UV absolute magnitudes
            
        Returns
        -------
        phi : ndarray
            Luminosity function Φ(M_UV) in Mpc^-3 mag^-1
        """
        return compute_luminosity_function(
            MUV_arr, self.z,
            [self.eps0, np.log10(self.Mc), self.a, self.b],
            self.nsat_params,
            self.sigma_UV,
            self.add_dust
        )
    
    def compute_bias(self, MUV_arr):
        """Compute galaxy bias.
        
        Parameters
        ----------
        MUV_arr : array_like
            Array of UV absolute magnitudes
            
        Returns
        -------
        bias : ndarray
            Galaxy bias as function of M_UV
        """
        return compute_galaxy_bias(
            MUV_arr, self.z,
            [self.eps0, np.log10(self.Mc), self.a, self.b],
            self.nsat_params,
            self.sigma_UV,
            self.add_dust
        )
    
    def mean_halo_mass(self, MUV_thresh):
        """Compute mean halo mass for galaxies brighter than threshold.
        
        Parameters
        ----------
        MUV_thresh : float
            UV magnitude threshold
            
        Returns
        -------
        log10_Mh : float
            Mean log10 halo mass
        """
        return compute_mean_halo_mass(
            MUV_thresh, self.z,
            [self.eps0, np.log10(self.Mc), self.a, self.b],
            self.nsat_params,
            self.sigma_UV,
            self.add_dust
        )
    
    def mean_bias(self, MUV_thresh):
        """Compute mean galaxy bias for galaxies brighter than threshold.
        
        Parameters
        ----------
        MUV_thresh : float
            UV magnitude threshold
            
        Returns
        -------
        bias : float
            Mean galaxy bias
        """
        return compute_mean_bias(
            MUV_thresh, self.z,
            [self.eps0, np.log10(self.Mc), self.a, self.b],
            self.nsat_params,
            self.sigma_UV,
            self.add_dust
        )
    
    def update_parameters(self, **kwargs):
        """Update model parameters.
        
        Parameters
        ----------
        **kwargs : dict
            Parameter updates. Can include:
            - eps0, a, b : float
            - log_Mc : float (will be converted)
            - log_Mcut, log_Msat : float (will be converted)
            - asat : float
            - sigma_UV : float
            - z : float
            - add_dust : bool
        """
        if 'eps0' in kwargs:
            self.eps0 = kwargs['eps0']
        if 'log_Mc' in kwargs:
            self.Mc = 10**kwargs['log_Mc']
        elif 'Mc' in kwargs:
            self.Mc = kwargs['Mc']
        if 'a' in kwargs:
            self.a = kwargs['a']
        if 'b' in kwargs:
            self.b = kwargs['b']
        if 'log_Mcut' in kwargs:
            self.Mcut = 10**kwargs['log_Mcut']
        elif 'Mcut' in kwargs:
            self.Mcut = kwargs['Mcut']
        if 'log_Msat' in kwargs:
            self.Msat = 10**kwargs['log_Msat']
        elif 'Msat' in kwargs:
            self.Msat = kwargs['Msat']
        if 'asat' in kwargs:
            self.asat = kwargs['asat']
        if 'sigma_UV' in kwargs:
            self.sigma_UV = kwargs['sigma_UV']
        if 'z' in kwargs:
            self.z = kwargs['z']
            # Clear cached HMF since redshift changed
            self._hmf = None
            self._log10_mass = None
        if 'add_dust' in kwargs:
            self.add_dust = kwargs['add_dust']
    
    def __repr__(self):
        return (f"UVLuminosityFunction(z={self.z}, "
                f"eps0={self.eps0:.3f}, Mc={self.Mc:.2e}, "
                f"a={self.a:.2f}, b={self.b:.2f})")
    
    def __str__(self):
        s = f"UV Luminosity Function Model at z={self.z}\n"
        s += "="*50 + "\n"
        s += "UVHMR Parameters:\n"
        s += f"  eps0 = {self.eps0:.3f}\n"
        s += f"  Mc = {self.Mc:.2e} M_sun\n"
        s += f"  a = {self.a:.2f}\n"
        s += f"  b = {self.b:.2f}\n"
        s += "\nSatellite Parameters:\n"
        s += f"  Mcut = {self.Mcut:.2e} M_sun\n"
        s += f"  Msat = {self.Msat:.2e} M_sun\n"
        s += f"  asat = {self.asat:.2f}\n"
        s += "\nOther:\n"
        s += f"  sigma_UV = {self.sigma_UV:.2f} mag\n"
        s += f"  add_dust = {self.add_dust}\n"
        s += "="*50
        return s