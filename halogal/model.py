"""Unified model classes for UVHMR and HOD calculations.

This module provides a class-based interface for all galaxy population
modeling calculations, combining UV-halo mass relations with halo
occupation distributions.
"""

import numpy as np
from scipy import interpolate, integrate
from scipy.special import erf
from scipy.stats import norm

from .config import BARYON_FRACTION
from .luminosity import MUV_from_SFR
from .cosmology import get_halo_mass_function, get_halo_bias

try:
    import halomod
    import halomod.hod
    HALOMOD_AVAILABLE = True
except ImportError:
    HALOMOD_AVAILABLE = False


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
    galaxy occupation functions.
    
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
    >>> # Use all defaults
    >>> model = HODModel(z=6.0)
    >>> 
    >>> # Compute occupation numbers
    >>> import numpy as np
    >>> Mh = np.logspace(10, 13, 100)
    >>> N_cen = model.Ncen(Mh, MUV_thresh=-19.0)
    >>> N_sat = model.Nsat(Mh, MUV_thresh=-19.0)
    """
    
    def __init__(self, z, eps0=None, Mc=None, a=None, b=None, 
                 sigma_UV=None, Mcut=None, Msat=None, asat=None, 
                 add_dust=True):
        """Initialize HOD model."""
        from .config import DEFAULT_HOD_PARAMS
        
        # Initialize parent class
        super().__init__(z, eps0, Mc, a, b, add_dust)
        
        # Set HOD parameters
        self.sigma_UV = sigma_UV if sigma_UV is not None else DEFAULT_HOD_PARAMS['sigma_UV']
        self.Mcut = Mcut if Mcut is not None else 10**DEFAULT_HOD_PARAMS['Mcut']
        self.Msat = Msat if Msat is not None else 10**DEFAULT_HOD_PARAMS['Msat']
        self.asat = asat if asat is not None else DEFAULT_HOD_PARAMS['asat']
    
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



class Observables:
    """Compute observables from an HOD model.
    
    This class provides methods to compute various observables like
    luminosity functions, galaxy bias, and correlation functions from
    an HOD model.
    
    All observable methods accept optional ``**params`` keyword arguments
    (eps0, Mc, a, b, sigma_UV, Mcut, Msat, asat) to update the underlying
    HOD model before computing. This enables efficient MCMC loops without
    recreating any objects.
    
    For correlation functions specifically, use the initialize/update pattern:
    1. Call initialize_correlation_model() once
    2. Call update_correlation_model() in your MCMC loop
    
    Parameters
    ----------
    hod_model : HODModel
        The HOD model to use for computations
        
    Attributes
    ----------
    hod_model : HODModel
        The underlying HOD model
        
    Examples
    --------
    >>> # Basic usage
    >>> model = HODModel(z=6.0)
    >>> obs = Observables(model)
    >>> MUV = np.linspace(-22, -16, 20)
    >>> phi = obs.luminosity_function(MUV)
    
    >>> # MCMC fitting - all observables support inline parameter updates
    >>> obs = Observables(HODModel(z=6.0))
    >>> for eps0, sigma_UV in mcmc_samples:
    ...     phi = obs.luminosity_function(MUV, eps0=eps0, sigma_UV=sigma_UV)
    ...     bg = obs.galaxy_bias(MUV, eps0=eps0, sigma_UV=sigma_UV)
    
    >>> # Correlation function MCMC (uses halomod caching)
    >>> obs.initialize_correlation_model(MUV_thresh1=-19.1)
    >>> for eps0, sigma_UV in mcmc_samples:
    ...     result = obs.update_correlation_model(eps0=eps0, sigma_UV=sigma_UV)
    ...     w_theta = result['correlation']
    """
    
    # HOD parameter names accepted by update_parameters
    _HOD_PARAM_NAMES = frozenset([
    'z', 'eps0', 'Mc', 'a', 'b', 'sigma_UV', 'Mcut', 'Msat', 'asat', 'add_dust'
])
    
    def __init__(self, hod_model):
        """Initialize observables calculator.
        
        Parameters
        ----------
        hod_model : HODModel
            The HOD model to use
        """
        self.hod_model = hod_model
        self._hmf_cache = {}
        
        # For efficient correlation function updates
        self._halomod_model = None
        self._halomod_config = {}
        self._hod_wrapper = None
    
    def _apply_params(self, params):
        """Apply parameter updates to the underlying HOD model.
        
        Extracts recognized HOD parameters from ``params``, applies them
        to ``self.hod_model``, and returns any remaining kwargs.
        
        Parameters
        ----------
        params : dict
            Keyword arguments, possibly containing HOD parameter overrides.
            
        Returns
        -------
        remaining : dict
            Any kwargs that are not HOD parameters (passed through).
        """
        hod_updates = {k: v for k, v in params.items() if k in self._HOD_PARAM_NAMES}
        remaining = {k: v for k, v in params.items() if k not in self._HOD_PARAM_NAMES}
        if hod_updates:
            self.hod_model.update_parameters(**hod_updates)
        return remaining

    def _get_hmf(self, recalculate=False):
        """Get HMF - from halomod if available, otherwise from Colossus.
        Use halomod HMF when available for
        consistency with correlation functions, otherwise use Colossus.
        Parameters
        ----------
        recalculate : bool, optional
            Force recalculation. Default is False.
        Returns
        -------
        log10_mass : ndarray
            Array of log10(M/M_sun) in physical units
        hmf : ndarray
            Halo mass function dn/dln(M) in physical units of Mpc^-3
        """
        # If halomod model exists, use its HMF for consistency
        if self._halomod_model is not None:
            # Extract from halomod
            M_h = self._halomod_model.m  # M_sun/h
            dndlnm_h = self._halomod_model.dndlnm  # (h/Mpc)^3

            # Convert to physical units
            h = self._halomod_config['h']
            M_physical = M_h / h  # M_sun
            hmf_physical = dndlnm_h * h**3  # Mpc^-3

            return np.log10(M_physical), hmf_physical

        # Otherwise use Colossus (backward compatible)
        cache_key = self.hod_model.z
        if recalculate or cache_key not in self._hmf_cache:
            log10_Mh, hmf = get_halo_mass_function(
                self.hod_model.z, M_min=6, M_max=15, num_points=2048
            )
            self._hmf_cache[cache_key] = (log10_Mh, hmf)
        return self._hmf_cache[cache_key]
    
    
    def luminosity_function(self, MUV_arr, dMUV=0.05, **params):
        """Compute UV luminosity function.
        
        Parameters
        ----------
        MUV_arr : array_like
            Array of UV absolute magnitudes
        dMUV : float, optional
            Magnitude bin width. Default is 0.05.
        **params : dict, optional
            HOD parameter overrides (eps0, Mc, a, b, sigma_UV, Mcut,
            Msat, asat, add_dust). Applied to the model before computing.
            
        Returns
        -------
        phi : ndarray
            Luminosity function Φ(M_UV) in Mpc^-3 mag^-1
            
        Examples
        --------
        >>> obs.luminosity_function(MUV)                          # current params
        >>> obs.luminosity_function(MUV, eps0=0.3, sigma_UV=0.5)  # override
        """
        self._apply_params(params)
        
        MUV_arr = np.atleast_1d(MUV_arr)
        log10_Mh, hmf = self._get_hmf()
        Mh = 10**log10_Mh
        
        phi = np.zeros_like(MUV_arr)
        
        for i, muv in enumerate(MUV_arr):
            # Number in magnitude bin
            N_cen_bin = (self.hod_model.Ncen(Mh, muv + dMUV/2) - 
                        self.hod_model.Ncen(Mh, muv - dMUV/2))
            N_sat_bin = (self.hod_model.Nsat(Mh, muv + dMUV/2) - 
                        self.hod_model.Nsat(Mh, muv - dMUV/2))
            
            integrand = (N_cen_bin + N_sat_bin) * hmf
            phi[i] = integrate.simpson(integrand, x=log10_Mh) / dMUV
        
        return phi
    
    def galaxy_bias(self, MUV_arr, dMUV=0.05, **params):
        """Compute galaxy bias as a function of UV magnitude.
        
        Parameters
        ----------
        MUV_arr : array_like
            Array of UV absolute magnitudes
        dMUV : float, optional
            Magnitude bin width
        **params : dict, optional
            HOD parameter overrides (eps0, Mc, a, b, sigma_UV, Mcut,
            Msat, asat, add_dust). Applied to the model before computing.
            
        Returns
        -------
        bias : ndarray
            Galaxy bias
        """
        self._apply_params(params)
        
        MUV_arr = np.atleast_1d(MUV_arr)
        log10_Mh, hmf = self._get_hmf()
        Mh = 10**log10_Mh
        bias_halo = get_halo_bias(Mh, self.hod_model.z, mdef='vir', model='tinker10')
        
        bg = np.zeros_like(MUV_arr)
        
        for i, muv in enumerate(MUV_arr):
            N_cen_bin = (self.hod_model.Ncen(Mh, muv + dMUV/2) - 
                        self.hod_model.Ncen(Mh, muv - dMUV/2))
            N_sat_bin = (self.hod_model.Nsat(Mh, muv + dMUV/2) - 
                        self.hod_model.Nsat(Mh, muv - dMUV/2))
            
            bg_integrand = (N_cen_bin + N_sat_bin) * hmf * bias_halo
            ngal_integrand = (N_cen_bin + N_sat_bin) * hmf
            
            ngal = integrate.simpson(ngal_integrand, x=log10_Mh)
            
            if ngal > 0:
                bg[i] = integrate.simpson(bg_integrand, x=log10_Mh) / ngal
            else:
                bg[i] = 0.0
        
        return bg
    
    def mean_halo_mass(self, MUV_thresh, **params):
        """Compute mean halo mass for galaxies brighter than threshold.
        
        Parameters
        ----------
        MUV_thresh : float
            UV magnitude threshold
        **params : dict, optional
            HOD parameter overrides.
            
        Returns
        -------
        log10_Mh_mean : float
            Mean log10 halo mass
        """
        self._apply_params(params)
        
        log10_Mh, hmf = self._get_hmf()
        Mh = 10**log10_Mh
        
        N_gal = self.hod_model.Ngal(Mh, MUV_thresh)
        
        mh_integrand = N_gal * hmf * Mh
        ngal_integrand = N_gal * hmf
        
        ngal = integrate.simpson(ngal_integrand, x=log10_Mh)
        
        if ngal > 0:
            mh_mean = integrate.simpson(mh_integrand, x=log10_Mh) / ngal
            return np.log10(mh_mean)
        else:
            return np.nan
    
    def mean_bias(self, MUV_thresh, **params):
        """Compute mean galaxy bias for galaxies brighter than threshold.
        
        Parameters
        ----------
        MUV_thresh : float
            UV magnitude threshold
        **params : dict, optional
            HOD parameter overrides.
            
        Returns
        -------
        bias_mean : float
            Mean galaxy bias
        """
        self._apply_params(params)
        
        log10_Mh, hmf = self._get_hmf()
        Mh = 10**log10_Mh
        bias_halo = get_halo_bias(Mh, self.hod_model.z, mdef='vir', model='tinker10')
        
        N_gal = self.hod_model.Ngal(Mh, MUV_thresh)
        
        bg_integrand = N_gal * hmf * bias_halo
        ngal_integrand = N_gal * hmf
        
        ngal = integrate.simpson(ngal_integrand, x=log10_Mh)
        
        if ngal > 0:
            return integrate.simpson(bg_integrand, x=log10_Mh) / ngal
        else:
            return np.nan
    
    def number_density(self, MUV_thresh, **params):
        """Compute cumulative galaxy number density brighter than threshold.
        
        Parameters
        ----------
        MUV_thresh : float
            UV magnitude threshold
        **params : dict, optional
            HOD parameter overrides.
            
        Returns
        -------
        n_gal : float
            Galaxy number density in Mpc^-3
        """
        self._apply_params(params)
        
        log10_Mh, hmf = self._get_hmf()
        Mh = 10**log10_Mh
        
        N_gal = self.hod_model.Ngal(Mh, MUV_thresh)
        integrand = N_gal * hmf
        
        return integrate.simpson(integrand, x=log10_Mh)
    
    def initialize_correlation_model(self, MUV_thresh1, MUV_thresh2=0,
                                    cosmo_model=None, correlation_type='angular',
                                    zmin=None, zmax=None, znum=150,
                                    p1=None,
                                    theta_min=1.0, theta_max=7200.0, theta_num=100,
                                    rmin=1e-3, rmax=200.0, rnum=200,
                                    rp_min=0.1, rp_max=100.0, rp_num=100,
                                    pi_max=100.0,
                                    hmf_model='Watson', bias_model='Tinker10',
                                    concentration_model='Duffy08',
                                    halo_profile='NFW',
                                    exclusion_model='DblSphere_',
                                    sd_bias_model=None,
                                    **kwargs):
        """Initialize halomod model for efficient parameter updates.
        
        Call this once before MCMC fitting. Then use update_correlation_model()
        to efficiently update HOD parameters.
        
        This creates and caches the halomod model, allowing fast parameter
        updates via halomod's update() method.
        
        Parameters are the same as compute_correlation_function().
        
        Returns
        -------
        result : dict
            Dictionary with 'model', 'separation', 'correlation', 'type'
            
        Examples
        --------
        >>> # Initialize once
        >>> obs = Observables(model)
        >>> result = obs.initialize_correlation_model(
        ...     MUV_thresh1=-19.1, correlation_type='angular'
        ... )
        >>> 
        >>> # Update efficiently in MCMC loop
        >>> for eps0 in mcmc_chain:
        ...     result = obs.update_correlation_model(eps0=eps0)
        ...     chi2 = compute_chi2(result['correlation'], data)
        """
        if not HALOMOD_AVAILABLE:
            raise ImportError("halomod package required")
        
        # Set up redshift range
        z_center = self.hod_model.z
        if zmin is None:
            zmin = max(0.1, z_center - 1.0)
        if zmax is None:
            zmax = z_center + 1.0
        
        # Set up cosmology and extract h
        if cosmo_model is None:
            from .config import DEFAULT_COSMO_PARAMS
            cosmo_params = {
                'Om0': DEFAULT_COSMO_PARAMS['Om0'],
                'Ob0': DEFAULT_COSMO_PARAMS['Ob0'],
                'H0': DEFAULT_COSMO_PARAMS['H0'],
            }
            h = DEFAULT_COSMO_PARAMS['H0'] / 100.0
        elif isinstance(cosmo_model, dict):
            cosmo_params = cosmo_model
            if 'H0' not in cosmo_params:
                raise ValueError("cosmo_model dict must include 'H0'")
            h = cosmo_params['H0'] / 100.0
        else:
            cosmo_params = cosmo_model
            h = cosmo_model.H0.value / 100.0
        
        # Store configuration for updates
        self._halomod_config = {
            'MUV_thresh1': MUV_thresh1,
            'MUV_thresh2': MUV_thresh2,
            'correlation_type': correlation_type,
            'h': h,
            'z_center': z_center,
            'zmin': zmin,
            'zmax': zmax,
            'znum': znum,
            'cosmo_params': cosmo_params,
        }


        # Halomod parameters
        halomod_params = {
            'cosmo_params' : cosmo_params,
            'z': z_center,
            'rmin': rmin * h,
            'rmax': rmax * h,
            'rnum': rnum,
            'Mmin': kwargs.pop('Mmin', 5),
            'Mmax': kwargs.pop('Mmax', 16),
            'dlog10m': kwargs.pop('dlog10m', 0.1),
            'dlnk': kwargs.pop('dlnk', 0.2),
            'lnk_min': kwargs.pop('lnk_min', -6),
            'lnk_max': kwargs.pop('lnk_max', 7),
            'hmf_model': hmf_model,
            'tracer_concentration_model': concentration_model,
            'halo_profile_model': halo_profile,
            'hc_spectrum': kwargs.pop('hc_spectrum', 'nonlinear'),
            'exclusion_model': exclusion_model,
            'bias_model': bias_model,
            'sd_bias_model': sd_bias_model,
            'hod_model' : HalomodHOD,
            'hod_params' : {
                'z' : z_center,
                'MUV_thresh1' : MUV_thresh1,
                'MUV_thresh2' : MUV_thresh2,
                'eps0' : self.hod_model.eps0,
                'Mc' : self.hod_model.Mc,
                'a' : self.hod_model.a,
                'b' : self.hod_model.b,
                'sigma_UV' : self.hod_model.sigma_UV,
                'Mcut' : self.hod_model.Mcut,
                'Msat' : self.hod_model.Msat,
                'asat' : self.hod_model.asat,
                'add_dust' : self.hod_model.add_dust,
                'h' : h,
                'M_min': kwargs.pop('M_min', 9.0)}
        }

        halomod_params.update(kwargs)


        if correlation_type.lower() == 'angular':
            theta_min_rad = theta_min / 3600.0 * (np.pi / 180.0)
            theta_max_rad = theta_max / 3600.0 * (np.pi / 180.0)
            
            if p1 is None:
                # Create a normalized callable directly
                def p1(z):
                    return norm.pdf(z, loc=z_center, scale=0.5)
            
            self._halomod_model = halomod.AngularCF(
                zmin=zmin,
                zmax=zmax,
                znum=znum,
                p1=p1,
                theta_min=theta_min_rad,
                theta_max=theta_max_rad,
                theta_num=theta_num,
                logu_min=kwargs.pop('logu_min', -3),
                logu_max=kwargs.pop('logu_max', 2.0),
                unum=kwargs.pop('unum', 100),
                dr_table=kwargs.pop('dr_table', 0.03),
                **halomod_params)
            
            separation = self._halomod_model.theta * (180.0 / np.pi) * 3600.0
            correlation = self._halomod_model.angular_corr_gal
            
        elif correlation_type.lower() == 'real':
            self._halomod_model = halomod.TracerHaloModel(
                **halomod_params
            )
            
            separation = self._halomod_model.r / h
            correlation = self._halomod_model.corr_auto_tracer
            
        elif correlation_type.lower() == 'projected':
            self._halomod_model = halomod.ProjectedCF(
                rp_min=rp_min * h,
                rp_max=rp_max * h,
                rp_num=rp_num,
                proj_limit=pi_max * h,
                **halomod_params
            )
            
            separation = self._halomod_model.r / h
            correlation = self._halomod_model.projected_corr_gal
            
        else:
            raise ValueError(f"Unknown correlation_type: {correlation_type}")
        
        return {
            'model': self._halomod_model,
            'separation': separation,
            'correlation': correlation,
            'type': correlation_type,
        }
    
    def update_correlation_model(self, z_center=None, p_of_z=None, MUV_thresh1=None, MUV_thresh2=None,
                                eps0=None, Mc=None, a=None, b=None,
                                sigma_UV=None, Mcut=None, Msat=None, asat=None):
        """Update HOD parameters and recompute correlation function efficiently.
        
        This uses halomod's update() method to efficiently recompute the
        correlation function with new parameters, without recreating the
        entire model.
        
        Must call initialize_correlation_model() first.
        
        Parameters
        ----------
        z_center : float, optional
            Update the redshift center
        p_of_z : callable, optional
            Update the redshift distribution function
        MUV_thresh1 : float, optional
            Update UV magnitude threshold
        MUV_thresh2 : float, optional
            Update faint-end threshold
        eps0 : float, optional
            Update star formation efficiency
        Mc : float, optional
            Update characteristic mass
        a : float, optional
            Update low-mass slope
        b : float, optional
            Update high-mass slope
        sigma_UV : float, optional
            Update UV scatter
        Mcut : float, optional
            Update satellite cutoff mass
        Msat : float, optional
            Update satellite normalization mass
        asat : float, optional
            Update satellite slope
            
        Returns
        -------
        result : dict
            Dictionary with 'model', 'separation', 'correlation', 'type'
            
        Examples
        --------
        >>> # Initialize once
        >>> result = obs.initialize_correlation_model(MUV_thresh1=-19.1)
        >>> 
        >>> # Update in MCMC loop
        >>> for eps0, sigma_UV in zip(eps0_chain, sigma_UV_chain):
        ...     result = obs.update_correlation_model(
        ...         eps0=eps0, sigma_UV=sigma_UV
        ...     )
        ...     likelihood = compute_likelihood(result['correlation'])
        """
        if self._halomod_model is None:
            raise ValueError(
                "Must call initialize_correlation_model() before update_correlation_model(). "
                "Or use compute_correlation_function() for one-off calculations."
            )
        
        # Update HODModel parameters
        update_dict = {}
        if eps0 is not None:
            update_dict['eps0'] = eps0
        if Mc is not None:
            update_dict['Mc'] = Mc
        if a is not None:
            update_dict['a'] = a
        if b is not None:
            update_dict['b'] = b
        if sigma_UV is not None:
            update_dict['sigma_UV'] = sigma_UV
        if Mcut is not None:
            update_dict['Mcut'] = Mcut
        if Msat is not None:
            update_dict['Msat'] = Msat
        if asat is not None:
            update_dict['asat'] = asat
        
        if update_dict:
            self.hod_model.update_parameters(**update_dict)
        
        # Update magnitude thresholds if provided
        if MUV_thresh1 is not None:
            self._halomod_config['MUV_thresh1'] = MUV_thresh1
        if MUV_thresh2 is not None:
            self._halomod_config['MUV_thresh2'] = MUV_thresh2
        
        # Update z_center if provided
        if z_center is not None:
            self._halomod_config['z_center'] = z_center
        
        # Build HOD params update dict (always include current values)
        hod_params_update = {
            'z': self._halomod_config['z_center'],
            'MUV_thresh1': self._halomod_config['MUV_thresh1'],
            'MUV_thresh2': self._halomod_config['MUV_thresh2'],
            'eps0': self.hod_model.eps0,
            'Mc': self.hod_model.Mc,
            'a': self.hod_model.a,
            'b': self.hod_model.b,
            'sigma_UV': self.hod_model.sigma_UV,
            'Mcut': self.hod_model.Mcut,
            'Msat': self.hod_model.Msat,
            'asat': self.hod_model.asat,
            'add_dust': self.hod_model.add_dust,
        }
        
        # Build halomod update dict (only include non-None values)
        halomod_update = {'hod_params': hod_params_update}
        
        if z_center is not None:
            halomod_update['z'] = z_center
        
        if p_of_z is not None:
            halomod_update['p1'] = p_of_z
        
        # Update halomod model efficiently
        self._halomod_model.update(**halomod_update)
        
        # Extract updated results
        correlation_type = self._halomod_config['correlation_type']
        h = self._halomod_config['h']
        
        if correlation_type.lower() == 'angular':
            separation = self._halomod_model.theta * (180.0 / np.pi) * 3600.0
            correlation = self._halomod_model.angular_corr_gal
        elif correlation_type.lower() == 'real':
            separation = self._halomod_model.r / h
            correlation = self._halomod_model.corr_auto_tracer
        elif correlation_type.lower() == 'projected':
            separation = self._halomod_model.rp / h
            correlation = self._halomod_model.projected_corr_gal
        
        return {
            'model': self._halomod_model,
            'separation': separation,
            'correlation': correlation,
            'type': correlation_type,
        }

    
    def compute_correlation_function(self, MUV_thresh1, MUV_thresh2=0, **kwargs):
        """Compute correlation function (convenience method).
        
        This creates a new halomod model each time. For MCMC fitting,
        use initialize_correlation_model() + update_correlation_model() instead.
        
        Parameters are the same as initialize_correlation_model().
        
        Returns
        -------
        result : dict
            Dictionary with 'model', 'separation', 'correlation', 'type'
        """
        return self.initialize_correlation_model(
            MUV_thresh1, MUV_thresh2, **kwargs
        )
    


class HalomodHOD(halomod.hod.HODPoisson if HALOMOD_AVAILABLE else object):
    """Halomod-compatible HOD class wrapping the HODModel.

    This class bridges our HODModel with halomod's expectations.
    
    IMPORTANT: Handles unit conversions between halomod (M_sun/h) and 
    HODModel (M_sun).

    Parameters
    ----------
    z : float
        Redshift
    MUV_thresh1 : float, optional
        Bright-end UV magnitude threshold
    MUV_thresh2 : float, optional
        Faint-end UV magnitude threshold
    eps0, Mc, a, b, sigma_UV, Mcut, Msat, asat, add_dust : optional
        HOD parameters. Default from config.
    h : float, optional
        Hubble parameter H0/100
        
    Examples
    --------
    >>> hod = HalomodHOD(z=5.0, MUV_thresh1=-19.1, h=0.7)
    """

    _defaults = {
        "z": 5.0,
        "MUV_thresh1": -16.0,
        "MUV_thresh2": 0.0,
        "eps0": None,
        "Mc": None,
        "a": None,
        "b": None,
        "sigma_UV": None,
        "Mcut": None,
        "Msat": None,
        "asat": None,
        "add_dust": True,
        "M_min": 8.0, 
    }

    def __init__(self, z=None, MUV_thresh1=None, MUV_thresh2=None,
                 eps0=None, Mc=None, a=None, b=None,
                 sigma_UV=None, Mcut=None, Msat=None, asat=None,
                 add_dust=None, h=None, M_min=None, **kwargs):
        """Initialize halomod HOD wrapper."""
        if not HALOMOD_AVAILABLE:
            raise ImportError("halomod package required")
        
        super().__init__(**kwargs)
        
        from .config import DEFAULT_HOD_PARAMS, DEFAULT_COSMO_PARAMS
        
        # Set parameters
        self.z = z if z is not None else self._defaults['z']
        self.MUV_thresh1 = MUV_thresh1 if MUV_thresh1 is not None else self._defaults['MUV_thresh1']
        self.MUV_thresh2 = MUV_thresh2 if MUV_thresh2 is not None else self._defaults['MUV_thresh2']
        
        self.eps0 = eps0 if eps0 is not None else DEFAULT_HOD_PARAMS['eps0']
        self.Mc = Mc if Mc is not None else 10**DEFAULT_HOD_PARAMS['logMc']
        self.a = a if a is not None else DEFAULT_HOD_PARAMS['a']
        self.b = b if b is not None else DEFAULT_HOD_PARAMS['b']
        self.sigma_UV = sigma_UV if sigma_UV is not None else DEFAULT_HOD_PARAMS['sigma_UV']
        self.Mcut = Mcut if Mcut is not None else 10**DEFAULT_HOD_PARAMS['Mcut']
        self.Msat = Msat if Msat is not None else 10**DEFAULT_HOD_PARAMS['Msat']
        self.asat = asat if asat is not None else DEFAULT_HOD_PARAMS['asat']
        self.add_dust = add_dust if add_dust is not None else DEFAULT_HOD_PARAMS['add_dust']
        self.M_min = M_min if M_min is not None else self._defaults['M_min']  # Add this
        
        self.h = h if h is not None else DEFAULT_COSMO_PARAMS['H0'] / 100.0
        
        # Create internal HODModel
        self.hod_model = HODModel(
            z=self.z, eps0=self.eps0, Mc=self.Mc, a=self.a, b=self.b,
            sigma_UV=self.sigma_UV, Mcut=self.Mcut, Msat=self.Msat,
            asat=self.asat, add_dust=self.add_dust
        )

    def _central_occupation(self, M):
        """Central occupation function for halomod.
        
        UNIT CONVERSION: M [M_sun] = M [M_sun/h] / h
        """
        M_physical = M / self.h
        
        if self.MUV_thresh2 > 0:
            N_c_bright = self.hod_model.Ncen(M_physical, self.MUV_thresh1)
            N_c_faint = self.hod_model.Ncen(M_physical, self.MUV_thresh2)
            return N_c_bright - N_c_faint
        else:
            return self.hod_model.Ncen(M_physical, self.MUV_thresh1)

    def _satellite_occupation(self, M):
        """Satellite occupation function for halomod.
        
        UNIT CONVERSION: M [M_sun] = M [M_sun/h] / h
        """
        M_physical = M / self.h
        
        if self.MUV_thresh2 > 0:
            N_s_bright = self.hod_model.Nsat(M_physical, self.MUV_thresh1)
            N_s_faint = self.hod_model.Nsat(M_physical, self.MUV_thresh2)
            return N_s_bright - N_s_faint
        else:
            return self.hod_model.Nsat(M_physical, self.MUV_thresh1)