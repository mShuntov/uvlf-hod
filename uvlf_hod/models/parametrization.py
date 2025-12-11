"""Redshift-dependent parametrization functions."""

import numpy as np


class RedshiftParametrization:
    """Class for handling redshift-dependent parameter evolution.
    
    Supports both linear and normalized redshift parametrizations:
    - Linear: param(z) = slope * z + offset
    - Normalized: param(z) = slope * z/(1+z) + offset
    """
    
    def __init__(self, use_linear=True):
        """Initialize parametrization scheme.
        
        Parameters
        ----------
        use_linear : bool, optional
            If True, use linear parametrization. If False, use z/(1+z).
            Default is True.
        """
        self.use_linear = use_linear
    
    def evaluate(self, z, slope, offset):
        """Evaluate parameter at given redshift.
        
        Parameters
        ----------
        z : float or array_like
            Redshift value(s)
        slope : float
            Evolution slope (d_param/dz or d_param/d[z/(1+z)])
        offset : float
            Value at z=0
            
        Returns
        -------
        param : float or array_like
            Parameter value at redshift z
        """
        z = np.atleast_1d(z)
        if self.use_linear:
            return slope * z + offset
        else:
            return slope * (z / (1 + z)) + offset


# Create default parametrization instance
default_param = RedshiftParametrization(use_linear=True)


def eps0_fz(z, deps_dz, eps_off, use_linear=True):
    """Star formation efficiency at redshift z.
    
    Parameters
    ----------
    z : float or array_like
        Redshift
    deps_dz : float
        Evolution slope
    eps_off : float
        Value at z=0
    use_linear : bool, optional
        Use linear vs normalized parametrization
        
    Returns
    -------
    eps0 : float or array_like
        Star formation efficiency parameter
    """
    param = RedshiftParametrization(use_linear=use_linear)
    return param.evaluate(z, deps_dz, eps_off)


def Mc_fz(z, dMc_dz, Mc_off, use_linear=True):
    """Characteristic halo mass at redshift z.
    
    Parameters
    ----------
    z : float or array_like
        Redshift
    dMc_dz : float
        Evolution slope
    Mc_off : float
        Value at z=0
    use_linear : bool, optional
        Use linear vs normalized parametrization
        
    Returns
    -------
    Mc : float or array_like
        Characteristic mass parameter
    """
    param = RedshiftParametrization(use_linear=use_linear)
    return param.evaluate(z, dMc_dz, Mc_off)


def a_fz(z, da_dz, a_off, use_linear=True):
    """Low-mass slope parameter at redshift z.
    
    Parameters
    ----------
    z : float or array_like
        Redshift
    da_dz : float
        Evolution slope
    a_off : float
        Value at z=0
    use_linear : bool, optional
        Use linear vs normalized parametrization
        
    Returns
    -------
    a : float or array_like
        Low-mass slope parameter
    """
    param = RedshiftParametrization(use_linear=use_linear)
    return param.evaluate(z, da_dz, a_off)


def b_fz(z, db_dz, b_off, use_linear=True):
    """High-mass slope parameter at redshift z.
    
    Parameters
    ----------
    z : float or array_like
        Redshift
    db_dz : float
        Evolution slope
    b_off : float
        Value at z=0
    use_linear : bool, optional
        Use linear vs normalized parametrization
        
    Returns
    -------
    b : float or array_like
        High-mass slope parameter
    """
    param = RedshiftParametrization(use_linear=use_linear)
    return param.evaluate(z, db_dz, b_off)


def sigma_UV_fz(z, dsigmaUV_dz, sigmaUV_off, use_linear=True):
    """UV magnitude scatter at redshift z.
    
    Parameters
    ----------
    z : float or array_like
        Redshift
    dsigmaUV_dz : float
        Evolution slope
    sigmaUV_off : float
        Value at z=0
    use_linear : bool, optional
        Use linear vs normalized parametrization
        
    Returns
    -------
    sigma_UV : float or array_like
        UV magnitude scatter
    """
    param = RedshiftParametrization(use_linear=use_linear)
    return param.evaluate(z, dsigmaUV_dz, sigmaUV_off)


def Mcut_fz(z, dMcut_dz, Mcut_off, use_linear=True):
    """Satellite cutoff mass at redshift z.
    
    Parameters
    ----------
    z : float or array_like
        Redshift
    dMcut_dz : float
        Evolution slope
    Mcut_off : float
        Value at z=0
    use_linear : bool, optional
        Use linear vs normalized parametrization
        
    Returns
    -------
    Mcut : float or array_like
        Cutoff mass parameter
    """
    param = RedshiftParametrization(use_linear=use_linear)
    return param.evaluate(z, dMcut_dz, Mcut_off)


def Msat_fz(z, dMsat_dz, Msat_off, use_linear=True):
    """Satellite normalization mass at redshift z.
    
    Parameters
    ----------
    z : float or array_like
        Redshift
    dMsat_dz : float
        Evolution slope
    Msat_off : float
        Value at z=0
    use_linear : bool, optional
        Use linear vs normalized parametrization
        
    Returns
    -------
    Msat : float or array_like
        Satellite mass parameter
    """
    param = RedshiftParametrization(use_linear=use_linear)
    return param.evaluate(z, dMsat_dz, Msat_off)


def asat_fz(z, dasat_dz, asat_off, use_linear=True):
    """Satellite power-law slope at redshift z.
    
    Parameters
    ----------
    z : float or array_like
        Redshift
    dasat_dz : float
        Evolution slope
    asat_off : float
        Value at z=0
    use_linear : bool, optional
        Use linear vs normalized parametrization
        
    Returns
    -------
    asat : float or array_like
        Satellite slope parameter
    """
    param = RedshiftParametrization(use_linear=use_linear)
    return param.evaluate(z, dasat_dz, asat_off)