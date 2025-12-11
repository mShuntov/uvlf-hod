"""Configuration and constants for uvlf-hod package."""

from astropy.cosmology import FlatLambdaCDM
from colossus.cosmology import cosmology

# Physical constants
BARYON_FRACTION = 0.16

# Default cosmological parameters
DEFAULT_COSMO_PARAMS = {
    'H0': 70.0,
    'Om0': 0.3,
    'Ob0': 0.0469,
    'Tcmb0': 2.725,
}

DEFAULT_COLOSSUS_PARAMS = {
    'flat': True,
    'H0': 70.0,
    'Om0': 0.3,
    'Ob0': 0.0469,
    'sigma8': 0.81,
    'ns': 1.0,
}

# Default HOD parameters
DEFAULT_HOD_PARAMS = {
    'z': 5.0,
    'MUV_thresh1': -16.0,
    'MUV_thresh2': 0.0,
    'eps0': 0.1,
    'logMc': 12.0,
    'alpha': 0.6,
    'beta': 0.35,
    'sigma_UV': 0.35,
    'Mcut': 8.0,
    'Msat': 12.0,
    'M_min': 9.0,
    'asat': 1.0,
    'add_dust': True,
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
    
    def __repr__(self):
        return f"CosmologyConfig(H0={self.params['H0']}, Om0={self.params['Om0']})"


# Default cosmology instance
default_cosmology = CosmologyConfig()
