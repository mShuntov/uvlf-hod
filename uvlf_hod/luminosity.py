"""UV luminosity and dust attenuation functions."""

import numpy as np

from .config import (
    BETA_Z_DATA, BETA_MUV_AT_M0, DBETA_DMUV, BETA_MUV0, BETA_C,
    DUST_C0, DUST_C1, DUST_SIGMA_BETA, DUST_Z_MAX,
    C_UV_DEFAULT
)


def beta_color(z, MUV):
    """UV spectral slope (beta) as a function of redshift and magnitude.
    
    Based on Bouwens+2013-14 data.
    
    Parameters
    ----------
    z : float or array_like
        Redshift
    MUV : float or array_like
        UV absolute magnitude
        
    Returns
    -------
    beta : float or ndarray
        UV spectral slope
        
    References
    ----------
    Bouwens et al. 2013, 2014
    """
    z = np.atleast_1d(z)
    MUV = np.atleast_1d(MUV)
    
    # Interpolate beta at M_UV = -19.5
    betaM0 = np.interp(z, BETA_Z_DATA, BETA_MUV_AT_M0, 
                       left=BETA_MUV_AT_M0[0], right=BETA_MUV_AT_M0[-1])
    
    # Interpolate slope
    dbetaM0 = (MUV - BETA_MUV0).T * np.interp(
        z, BETA_Z_DATA, DBETA_DMUV,
        left=DBETA_DMUV[0], right=DBETA_DMUV[-1]
    )
    
    # Two solutions depending on whether MUV > MUV0
    sol_bright = (betaM0 - BETA_C) * np.exp(dbetaM0 / (betaM0 - BETA_C)) + BETA_C
    sol_faint = dbetaM0 + betaM0
    
    # Apply Heaviside step function
    result = (sol_bright.T * np.heaviside(MUV - BETA_MUV0, 0.5) +
              sol_faint.T * np.heaviside(BETA_MUV0 - MUV, 0.5))
    
    return result


def dust_attenuation(z, MUV, high_z_dust=True, C0=None, C1=None):
    """UV dust attenuation as a function of redshift and magnitude.
    
    Parameters
    ----------
    z : float or array_like
        Redshift
    MUV : float or array_like
        UV absolute magnitude (observed)
    high_z_dust : bool, optional
        Whether to include dust at z > 8. Default is True.
    C0 : float, optional
        Attenuation constant. Default from config.
    C1 : float, optional
        Attenuation slope. Default from config.
        
    Returns
    -------
    A_UV : float or ndarray
        UV dust attenuation in magnitudes
        
    Notes
    -----
    This should be applied iteratively since the MUV-beta relation
    is calibrated on observations, not intrinsic magnitudes.
    """
    if C0 is None:
        C0 = DUST_C0
    if C1 is None:
        C1 = DUST_C1
    
    z = np.atleast_1d(z)
    MUV = np.atleast_1d(MUV)
    
    beta_z_muv = beta_color(z, MUV)
    
    # Calculate attenuation
    Auv = (C0 + 0.2 * np.log(10) * DUST_SIGMA_BETA**2 * C1**2 + 
           C1 * beta_z_muv)
    
    Auv = Auv.T
    
    # Turn off dust at high-z if requested
    if not high_z_dust:
        Auv *= np.heaviside(DUST_Z_MAX - z, 0.5)
    
    Auv = Auv.T
    
    # Ensure non-negative
    return np.maximum(Auv, 0.0)


def MUV_from_SFR(sfr, z, add_dust=True, c_uv=None, convergence_threshold=0.02):
    """Convert star formation rate to UV absolute magnitude.
    
    Parameters
    ----------
    sfr : float or array_like
        Star formation rate in M_sun/yr
    z : float
        Redshift
    add_dust : bool, optional
        Whether to add dust attenuation. Default is True.
    c_uv : float, optional
        SFR to UV luminosity conversion factor.
        Default is Salpeter IMF value from config.
    convergence_threshold : float, optional
        Convergence criterion for dust iteration. Default is 0.02 (2%).
        
    Returns
    -------
    MUV_obs : float or ndarray
        Observed UV absolute magnitude (with dust if add_dust=True)
        
    Notes
    -----
    Dust attenuation is added iteratively until convergence because
    the MUV-beta relation is calibrated on observations.
    """
    if c_uv is None:
        c_uv = C_UV_DEFAULT
    
    sfr = np.atleast_1d(sfr)
    
    # Convert SFR to UV luminosity
    LUV = sfr / c_uv  # erg/s/Hz
    
    # Intrinsic AB magnitude
    M_UV = 51.63 - 2.5 * np.log10(LUV)
    
    # Add dust attenuation iteratively
    M_UV_obs = M_UV.copy()
    
    if add_dust:
        M_UV_obs_prev = np.ones_like(M_UV_obs)
        max_iterations = 20
        iteration = 0
        
        while (np.sum(np.abs((M_UV_obs_prev - M_UV_obs) / M_UV_obs)) > 
               convergence_threshold and iteration < max_iterations):
            M_UV_obs_prev = M_UV_obs.copy()
            M_UV_obs = M_UV + dust_attenuation(z, M_UV_obs)
            iteration += 1
    
    return M_UV_obs


def SFR_from_MUV(MUV, z, remove_dust=True, c_uv=None):
    """Convert UV absolute magnitude to star formation rate.
    
    Parameters
    ----------
    MUV : float or array_like
        UV absolute magnitude
    z : float
        Redshift
    remove_dust : bool, optional
        Whether to remove dust attenuation. Default is True.
    c_uv : float, optional
        SFR to UV luminosity conversion factor.
        
    Returns
    -------
    sfr : float or ndarray
        Star formation rate in M_sun/yr
    """
    if c_uv is None:
        c_uv = C_UV_DEFAULT
    
    MUV = np.atleast_1d(MUV)
    
    # Remove dust if requested
    if remove_dust:
        MUV_intrinsic = MUV - dust_attenuation(z, MUV)
    else:
        MUV_intrinsic = MUV
    
    # Convert to luminosity
    LUV = 10**((51.63 - MUV_intrinsic) / 2.5)
    
    # Convert to SFR
    sfr = LUV * c_uv
    
    return sfr
