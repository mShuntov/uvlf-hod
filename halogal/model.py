"""Unified model classes for UVHMR and HOD calculations.

This module provides a class-based interface for all galaxy population
modeling calculations, combining UV-halo mass relations with halo
occupation distributions.
"""

import numpy as np
from scipy import interpolate, integrate
from scipy.special import erf

from .config import BARYON_FRACTION
from .luminosity import MUV_from_SFR
from .cosmology import get_halo_mass_function, get_halo_bias


class UVHMRModel:
    """Base class for UV-Halo Mass Relation models.
    
    This class handles the relationship between halo mass and UV luminosity
    through star formation. It serves as the foundation for more complex
    models that include occupation distributions.
    
    Parameters
    ----------
    z : float
        Redshift
    eps0 : float, optional
        Star formation efficiency normalization. Default from config.
    Mc : float, optional
        Characteristic halo mass in M_sun. Default from config.
    a : float, optional
        Low-mass power-law slope. Default from config.
    b : float, optional
        High-mass power-law slope. Default from config.
    add_dust : bool, optional
        Whether to include dust attenuation. Default is True.
        
    Attributes
    ----------
    z, eps0, Mc, a, b : float
        Model parameters
    add_dust : bool
        Whether dust is included
        
    Examples
    --------
    >>> # Use all defaults (fitted to Bouwens+2021)
    >>> model = UVHMRModel(z=6.0)
    >>> 
    >>> # Override specific parameters
    >>> model = UVHMRModel(z=6.0, eps0=0.2, Mc=10**12)
    """
    
    def __init__(self, z, eps0=None, Mc=None, a=None, b=None, add_dust=True):
        """Initialize UVHMR model."""
        from .config import DEFAULT_HOD_PARAMS
        
        self.z = z
        self.eps0 = eps0 if eps0 is not None else DEFAULT_HOD_PARAMS['eps0']
        self.Mc = Mc if Mc is not None else 10**DEFAULT_HOD_PARAMS['logMc']
        self.a = a if a is not None else DEFAULT_HOD_PARAMS['a']
        self.b = b if b is not None else DEFAULT_HOD_PARAMS['b']
        self.add_dust = add_dust
    
    def star_formation_efficiency(self, Mh):
        """Compute star formation efficiency.
        
        Parameters
        ----------
        Mh : float or array_like
            Halo mass in M_sun
            
        Returns
        -------
        epsilon : float or ndarray
            Star formation efficiency
        """
        Mh = np.atleast_1d(Mh)
        return 2 * self.eps0 / ((Mh / self.Mc)**(-self.a) + (Mh / self.Mc)**self.b)
    
    def halo_accretion_rate(self, Mh):
        """Compute halo mass accretion rate.
        
        Parameters
        ----------
        Mh : float or array_like
            Halo mass in M_sun
            
        Returns
        -------
        dMh_dt : float or ndarray
            Halo mass accretion rate in M_sun/yr
        """
        Mh = np.atleast_1d(Mh)
        return 0.03e-9 * Mh * (Mh / 1e12)**0.14 * (1 + self.z)**2.5
    
    def sfr(self, Mh):
        """Star formation rate from halo mass.
        
        Parameters
        ----------
        Mh : float or array_like
            Halo mass in M_sun
            
        Returns
        -------
        sfr : float or ndarray
            Star formation rate in M_sun/yr
        """
        f_star = self.star_formation_efficiency(Mh)
        return BARYON_FRACTION * f_star * self.halo_accretion_rate(Mh)
    
    def MUV(self, Mh):
        """UV absolute magnitude from halo mass.
        
        Parameters
        ----------
        Mh : float or array_like
            Halo mass in M_sun
            
        Returns
        -------
        MUV : float or ndarray
            UV absolute magnitude
        """
        sfr = self.sfr(Mh)
        return MUV_from_SFR(sfr, self.z, add_dust=self.add_dust)
    
    def Mhalo(self, MUV, M_min=1e6, M_max=1e16, num_points=1000):
        """Inverse: halo mass from UV magnitude.
        
        Parameters
        ----------
        MUV : float or array_like
            UV absolute magnitude
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
        Mhalo_grid = np.logspace(np.log10(M_min), np.log10(M_max), num_points)
        MUV_grid = self.MUV(Mhalo_grid)
        
        # Invert using interpolation
        interp_func = interpolate.interp1d(MUV_grid, Mhalo_grid, 
                                           fill_value='extrapolate')
        return interp_func(MUV)
    
    def update_parameters(self, **kwargs):
        """Update model parameters.
        
        Parameters
        ----------
        **kwargs : dict
            Parameters to update: z, eps0, Mc, a, b, add_dust
        """
        if 'z' in kwargs:
            self.z = kwargs['z']
        if 'eps0' in kwargs:
            self.eps0 = kwargs['eps0']
        if 'Mc' in kwargs:
            self.Mc = kwargs['Mc']
        if 'a' in kwargs:
            self.a = kwargs['a']
        if 'b' in kwargs:
            self.b = kwargs['b']
        if 'add_dust' in kwargs:
            self.add_dust = kwargs['add_dust']
    
    def __repr__(self):
        return (f"UVHMRModel(z={self.z}, eps0={self.eps0:.3f}, "
                f"Mc={self.Mc:.2e}, a={self.a:.2f}, b={self.b:.2f})")


class HODModel(UVHMRModel):
    """Halo Occupation Distribution model.
    
    This class extends UVHMRModel to include central and satellite
    galaxy occupation functions, enabling full galaxy population modeling.
    
    Parameters
    ----------
    z : float
        Redshift
    eps0 : float, optional
        Star formation efficiency normalization. Default from config.
    Mc : float, optional
        Characteristic halo mass in M_sun. Default from config.
    a : float, optional
        Low-mass power-law slope. Default from config.
    b : float, optional
        High-mass power-law slope. Default from config.
    sigma_UV : float, optional
        UV magnitude scatter. Default from config.
    Mcut : float, optional
        Satellite cutoff mass in M_sun. Default from config.
    Msat : float, optional
        Satellite normalization mass in M_sun. Default from config.
    asat : float, optional
        Satellite power-law slope. Default from config.
    add_dust : bool, optional
        Whether to include dust attenuation. Default is True.
        
    Attributes
    ----------
    sigma_UV, Mcut, Msat, asat : float
        HOD parameters
        
    Examples
    --------
    >>> # Use all defaults (fitted to Bouwens+2021 at z~5.4)
    >>> model = HODModel(z=6.0)
    >>> 
    >>> # Override specific parameters
    >>> model = HODModel(z=6.0, eps0=0.2, sigma_UV=0.5)
    >>> 
    >>> # Compute luminosity function
    >>> import numpy as np
    >>> MUV = np.linspace(-22, -16, 20)
    >>> phi = model.luminosity_function(MUV)
    """
    
    def __init__(self, z, eps0=None, Mc=None, a=None, b=None, 
                 sigma_UV=None, Mcut=None, Msat=None, asat=None, 
                 add_dust=True):
        """Initialize HOD model."""
        from .config import DEFAULT_HOD_PARAMS
        
        # Initialize parent class with defaults
        super().__init__(z, eps0, Mc, a, b, add_dust)
        
        # Set HOD parameters with defaults
        self.sigma_UV = sigma_UV if sigma_UV is not None else DEFAULT_HOD_PARAMS['sigma_UV']
        self.Mcut = Mcut if Mcut is not None else 10**DEFAULT_HOD_PARAMS['Mcut']
        self.Msat = Msat if Msat is not None else 10**DEFAULT_HOD_PARAMS['Msat']
        self.asat = asat if asat is not None else DEFAULT_HOD_PARAMS['asat']
        
        # Cache for HMF
        self._hmf_cache = {}
    
    def _get_hmf(self, recalculate=False):
        """Get or compute halo mass function."""
        cache_key = self.z
        if recalculate or cache_key not in self._hmf_cache:
            log10_Mh, hmf = get_halo_mass_function(
                self.z, M_min=6, M_max=15, num_points=2048
            )
            self._hmf_cache[cache_key] = (log10_Mh, hmf)
        return self._hmf_cache[cache_key]
    
    def Ncen(self, Mh, MUV_thresh):
        """Central galaxy occupation function.
        
        Average number of central galaxies with M_UV < MUV_thresh.
        
        Parameters
        ----------
        Mh : float or array_like
            Halo mass in M_sun
        MUV_thresh : float
            UV magnitude threshold (brighter than)
            
        Returns
        -------
        N_cen : float or ndarray
            Average number of central galaxies
        """
        MUV_mean = self.MUV(Mh)
        return 0.5 * (1 - erf((MUV_mean - MUV_thresh) / 
                              (np.sqrt(2) * self.sigma_UV)))
    
    def Nsat(self, Mh, MUV_thresh):
        """Satellite galaxy occupation function.
        
        Average number of satellite galaxies with M_UV < MUV_thresh.
        
        Parameters
        ----------
        Mh : float or array_like
            Halo mass in M_sun
        MUV_thresh : float
            UV magnitude threshold
            
        Returns
        -------
        N_sat : float or ndarray
            Average number of satellite galaxies
        """
        Mh = np.atleast_1d(Mh)
        N_cen = self.Ncen(Mh, MUV_thresh)
        
        # Satellites only exist above cutoff mass
        N_sat = np.zeros_like(Mh)
        mask = Mh > self.Mcut
        N_sat[mask] = ((Mh[mask] - self.Mcut) / self.Msat)**self.asat
        
        return N_sat * N_cen
    
    def Ngal(self, Mh, MUV_thresh):
        """Total galaxy occupation (central + satellite).
        
        Parameters
        ----------
        Mh : float or array_like
            Halo mass in M_sun
        MUV_thresh : float
            UV magnitude threshold
            
        Returns
        -------
        N_gal : float or ndarray
            Average total number of galaxies
        """
        return self.Ncen(Mh, MUV_thresh) + self.Nsat(Mh, MUV_thresh)
    
    def luminosity_function(self, MUV_arr, dMUV=0.05):
        """Compute UV luminosity function.
        
        Parameters
        ----------
        MUV_arr : array_like
            Array of UV absolute magnitudes
        dMUV : float, optional
            Magnitude bin width. Default is 0.05.
            
        Returns
        -------
        phi : ndarray
            Luminosity function Φ(M_UV) in Mpc^-3 mag^-1
        """
        MUV_arr = np.atleast_1d(MUV_arr)
        log10_Mh, hmf = self._get_hmf()
        Mh = 10**log10_Mh
        
        phi = np.zeros_like(MUV_arr)
        
        for i, muv in enumerate(MUV_arr):
            # Number in magnitude bin
            N_cen_bin = (self.Ncen(Mh, muv + dMUV/2) - 
                        self.Ncen(Mh, muv - dMUV/2))
            N_sat_bin = (self.Nsat(Mh, muv + dMUV/2) - 
                        self.Nsat(Mh, muv - dMUV/2))
            
            integrand = (N_cen_bin + N_sat_bin) * hmf
            phi[i] = integrate.simpson(integrand, x=log10_Mh) / dMUV
        
        return phi
    
    def galaxy_bias(self, MUV_arr, dMUV=0.05):
        """Compute galaxy bias as a function of UV magnitude.
        
        Parameters
        ----------
        MUV_arr : array_like
            Array of UV absolute magnitudes
        dMUV : float, optional
            Magnitude bin width
            
        Returns
        -------
        bias : ndarray
            Galaxy bias
        """
        MUV_arr = np.atleast_1d(MUV_arr)
        log10_Mh, hmf = self._get_hmf()
        Mh = 10**log10_Mh
        bias_halo = get_halo_bias(Mh, self.z, mdef='vir', model='tinker10')
        
        bg = np.zeros_like(MUV_arr)
        
        for i, muv in enumerate(MUV_arr):
            N_cen_bin = (self.Ncen(Mh, muv + dMUV/2) - 
                        self.Ncen(Mh, muv - dMUV/2))
            N_sat_bin = (self.Nsat(Mh, muv + dMUV/2) - 
                        self.Nsat(Mh, muv - dMUV/2))
            
            bg_integrand = (N_cen_bin + N_sat_bin) * hmf * bias_halo
            ngal_integrand = (N_cen_bin + N_sat_bin) * hmf
            
            ngal = integrate.simpson(ngal_integrand, x=log10_Mh)
            
            if ngal > 0:
                bg[i] = integrate.simpson(bg_integrand, x=log10_Mh) / ngal
            else:
                bg[i] = 0.0
        
        return bg
    
    def mean_halo_mass(self, MUV_thresh):
        """Compute mean halo mass for galaxies brighter than threshold.
        
        Parameters
        ----------
        MUV_thresh : float
            UV magnitude threshold
            
        Returns
        -------
        log10_Mh_mean : float
            Mean log10 halo mass
        """
        log10_Mh, hmf = self._get_hmf()
        Mh = 10**log10_Mh
        
        N_gal = self.Ngal(Mh, MUV_thresh)
        
        mh_integrand = N_gal * hmf * Mh
        ngal_integrand = N_gal * hmf
        
        ngal = integrate.simpson(ngal_integrand, x=log10_Mh)
        
        if ngal > 0:
            mh_mean = integrate.simpson(mh_integrand, x=log10_Mh) / ngal
            return np.log10(mh_mean)
        else:
            return np.nan
    
    def mean_bias(self, MUV_thresh):
        """Compute mean galaxy bias for galaxies brighter than threshold.
        
        Parameters
        ----------
        MUV_thresh : float
            UV magnitude threshold
            
        Returns
        -------
        bias_mean : float
            Mean galaxy bias
        """
        log10_Mh, hmf = self._get_hmf()
        Mh = 10**log10_Mh
        bias_halo = get_halo_bias(Mh, self.z, mdef='vir', model='tinker10')
        
        N_gal = self.Ngal(Mh, MUV_thresh)
        
        bg_integrand = N_gal * hmf * bias_halo
        ngal_integrand = N_gal * hmf
        
        ngal = integrate.simpson(ngal_integrand, x=log10_Mh)
        
        if ngal > 0:
            return integrate.simpson(bg_integrand, x=log10_Mh) / ngal
        else:
            return np.nan
    
    def update_parameters(self, **kwargs):
        """Update model parameters.
        
        Parameters
        ----------
        **kwargs : dict
            Parameters to update: z, eps0, Mc, a, b, sigma_UV, 
            Mcut, Msat, asat, add_dust
        """
        super().update_parameters(**kwargs)
        
        if 'sigma_UV' in kwargs:
            self.sigma_UV = kwargs['sigma_UV']
        if 'Mcut' in kwargs:
            self.Mcut = kwargs['Mcut']
        if 'Msat' in kwargs:
            self.Msat = kwargs['Msat']
        if 'asat' in kwargs:
            self.asat = kwargs['asat']
        
        # Clear HMF cache if redshift changed
        if 'z' in kwargs:
            self._hmf_cache.clear()
    
    def __repr__(self):
        return (f"HODModel(z={self.z}, eps0={self.eps0:.3f}, "
                f"Mc={self.Mc:.2e}, a={self.a:.2f}, b={self.b:.2f}, "
                f"sigma_UV={self.sigma_UV:.2f})")
    
    def __str__(self):
        s = f"HOD Model at z={self.z}\n"
        s += "="*50 + "\n"
        s += "UVHMR Parameters:\n"
        s += f"  eps0 = {self.eps0:.3f}\n"
        s += f"  Mc = {self.Mc:.2e} M_sun\n"
        s += f"  a = {self.a:.2f}\n"
        s += f"  b = {self.b:.2f}\n"
        s += "\nHOD Parameters:\n"
        s += f"  sigma_UV = {self.sigma_UV:.2f} mag\n"
        s += f"  Mcut = {self.Mcut:.2e} M_sun\n"
        s += f"  Msat = {self.Msat:.2e} M_sun\n"
        s += f"  asat = {self.asat:.2f}\n"
        s += "\nSettings:\n"
        s += f"  add_dust = {self.add_dust}\n"
        s += "="*50
        return s