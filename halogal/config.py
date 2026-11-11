"""Configuration and constants for uvlf-hod package."""

from astropy.cosmology import FlatLambdaCDM
from colossus.cosmology import cosmology

# Physical constants
BARYON_FRACTION = 0.16

# ---------------------------------------------------------------------------
# Cosmology: ONE source of truth
# ---------------------------------------------------------------------------
# Every consumer -- Colossus (halo mass function, halo bias), Astropy
# (distances) and halomod/hmf (correlation functions) -- is configured from
# this single dict. They used to be specified separately, which left the UVLF
# and the clustering describing different universes: Colossus ran with
# sigma8=0.81, ns=1.0 while halomod silently fell back to its own Planck18
# defaults of 0.8159 and 0.9667, a 30-48% difference in dn/dlnM at z=8.
#
# NOTE ON ns: 1.0 is not a measured value (Planck gives ~0.965) and looks like
# a placeholder. It is kept here because the fitted HOD parameters in
# DEFAULT_HOD_PARAMS were derived with it, so changing it silently would
# invalidate them. Change it in this one place when you are ready to refit.
DEFAULT_COSMO_PARAMS = {
    'H0': 70.0,
    'Om0': 0.3,
    'Ob0': 0.0469,
    'Tcmb0': 2.725,
    'sigma8': 0.81,
    'ns': 1.0,
}

# Colossus wants 'flat' and its own spelling; derived, never edited by hand.
DEFAULT_COLOSSUS_PARAMS = {
    'flat': True,
    'H0': DEFAULT_COSMO_PARAMS['H0'],
    'Om0': DEFAULT_COSMO_PARAMS['Om0'],
    'Ob0': DEFAULT_COSMO_PARAMS['Ob0'],
    'sigma8': DEFAULT_COSMO_PARAMS['sigma8'],
    'ns': DEFAULT_COSMO_PARAMS['ns'],
}


def halomod_cosmo_kwargs(params=None):
    """Translate halogal's cosmology into halomod/hmf keyword arguments.

    Parameters
    ----------
    params : dict, optional
        Cosmology parameters in halogal's spelling. Defaults to
        ``DEFAULT_COSMO_PARAMS``.

    Returns
    -------
    dict
        ``{'cosmo_params': {...}, 'sigma_8': ..., 'n': ..., 'transfer_model':
        ...}``, ready to splat into ``halomod.TracerHaloModel``.

    Notes
    -----
    halomod takes the background cosmology as ``cosmo_params`` (which override
    its Planck18 default) but the power-spectrum normalisation and tilt as the
    separate top-level arguments ``sigma_8`` and ``n``. Passing only the former
    is the trap this helper exists to close.

    ``transfer_model='EH'`` matches Colossus, whose default transfer function
    is Eisenstein & Hu (1998). halomod defaults to CAMB instead, and that
    mismatch alone left the two halo mass functions differing by 0.21 dex at
    z=8; with EH the residual is 0.077 dex, which is the remaining difference
    between Colossus's ``watson13`` and hmf's ``Watson`` implementations.
    """
    p = DEFAULT_COSMO_PARAMS if params is None else params
    return {
        'cosmo_params': {
            'Om0': p['Om0'],
            'Ob0': p['Ob0'],
            'H0': p['H0'],
        },
        'sigma_8': p['sigma8'],
        'n': p['ns'],
        'transfer_model': p.get('transfer_model', 'EH'),
    }

# Default HOD parameters (from fits at z~5.4)
# Based on Shuntov et al. 2025 2025A&A...699A.231S - see Table with posterior medians
DEFAULT_HOD_PARAMS = {
    'z': 5.4,
    'MUV_thresh1': -16.0,
    'MUV_thresh2': 0.0,
    'eps0': 0.19,        # Star formation efficiency
    'logMc': 11.64,      # log10(Mc/Msun) - characteristic mass
    'a': 0.69,           # Low-mass slope (beta in paper)
    'b': 0.65,           # High-mass slope (gamma in paper)
    'sigma_UV': 0.69,    # UV magnitude scatter
    'Mcut': 9.57,        # log10(Mcut/Msun) - satellite cutoff
    'Msat': 12.65,       # log10(Msat/Msun) - satellite normalization
    'M_min': 9.0,
    'asat': 0.85,        # Satellite power-law slope (alpha_sat in paper)
    'add_dust': True,
}

# Redshift evolution parameters (linear parametrization)
# param(z) = d_param_dz * z + C_param
DEFAULT_REDSHIFT_EVOLUTION = {
    'd_eps0_dz': 0.02,
    'C_eps0': 0.09,
    'd_logMc_dz': 0.14,
    'C_logMc': 10.86,
    'd_a_dz': 0.01,      # beta in paper
    'C_a': 0.62,
    'd_b_dz': -0.06,     # gamma in paper
    'C_b': 0.97,
    'd_sigmaUV_dz': -0.03,
    'C_sigmaUV': 0.74,
    'd_Mcut_dz': 0.84,
    'C_Mcut': 5.06,
    'd_Msat_dz': 0.06,
    'C_Msat': 12.38,
    'd_asat_dz': 0.04,
    'C_asat': 0.56,
}

# Beta-magnitude relation parameters (Bouwens 2013-14)
BETA_Z_DATA = [2.5, 3.8, 5.0, 5.9, 7.0, 8.0]
BETA_MUV_AT_M0 = [-1.7, -1.85, -1.91, -2.00, -2.05, -2.13]
DBETA_DMUV = [-0.20, -0.11, -0.14, -0.20, -0.20, -0.15]
BETA_MUV0 = -19.5
BETA_C = -2.33

# Dust attenuation parameters
DUST_C0 = 4.43
DUST_C1 = 1.99
DUST_SIGMA_BETA = 0.34
DUST_Z_MAX = 8.0

# SFR-to-UV conversion factor
C_UV_SALPETER = 1.15e-28  # (M_sun yr^-1) / (erg s^-1 Hz^-1)
C_UV_CHABRIER = 0.72e-28  # (M_sun yr^-1) / (erg s^-1 Hz^-1)

# Use Salpeter by default
C_UV_DEFAULT = C_UV_SALPETER

# Numerical parameters
HMF_NUM_POINTS = 2048
HMF_M_MIN = 6  # log10(M/M_sun)
HMF_M_MAX = 15  # log10(M/M_sun)


class CosmologyConfig:
    """Class to manage cosmology configurations."""
    
    def __init__(self, **kwargs):
        """Initialize cosmology with custom parameters."""
        self.params = {**DEFAULT_COSMO_PARAMS, **kwargs}
        self.colossus_params = {**DEFAULT_COLOSSUS_PARAMS, **kwargs}
        self._setup_cosmologies()
    
    def _setup_cosmologies(self):
        """Setup Astropy and Colossus cosmologies."""
        # Astropy cosmology
        self.astropy_cosmo = FlatLambdaCDM(
            H0=self.params['H0'],
            Om0=self.params['Om0'],
            Tcmb0=self.params.get('Tcmb0', 2.725),
            Ob0=self.params.get('Ob0', 0.0469)
        )
        
        # Colossus cosmology
        cosmology.setCosmology('custom', **self.colossus_params)
        self.colossus_cosmo = cosmology.getCurrent()

    def halomod_params(self):
        """Keyword arguments configuring halomod with THIS cosmology.

        Returns
        -------
        dict
            ``{'cosmo_params': {...}, 'sigma_8': ..., 'n': ...}``.

        Examples
        --------
        >>> cfg = CosmologyConfig(H0=67.4)
        >>> model = halomod.TracerHaloModel(z=8.0, **cfg.halomod_params())
        """
        return halomod_cosmo_kwargs(self.params)

    def __repr__(self):
        return (f"CosmologyConfig(H0={self.params['H0']}, "
                f"Om0={self.params['Om0']}, sigma8={self.params['sigma8']}, "
                f"ns={self.params['ns']})")


# Default cosmology instance
default_cosmology = CosmologyConfig()