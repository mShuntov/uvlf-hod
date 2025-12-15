"""Tests for core module."""

import pytest
import numpy as np
from halogal.core import (
    star_formation_efficiency,
    halo_accretion_rate,
    sfr_from_halo_mass,
    MUV_from_halo_mass,
    halo_mass_from_MUV,
    UVHaloMassRelation,
)


class TestStarFormationEfficiency:
    """Tests for star formation efficiency function."""
    
    def test_peak_at_Mc(self):
        """SFE should peak near Mc."""
        eps0, Mc, a, b = 0.1, 1e12, 0.6, 0.35
        
        # Evaluate at Mc
        sfe_Mc = star_formation_efficiency(Mc, eps0, Mc, a, b)
        
        # Should be close to 2*eps0 (exact value for double power law)
        assert np.isclose(sfe_Mc, 2 * eps0, rtol=0.01)
    
    def test_low_mass_powerlaw(self):
        """SFE should follow power law at low mass."""
        eps0, Mc, a, b = 0.1, 1e12, 0.6, 0.35
        
        # Test at low mass
        M_low = Mc / 100
        sfe_low = star_formation_efficiency(M_low, eps0, Mc, a, b)
        
        # Should be much smaller than peak
        assert sfe_low < eps0
    
    def test_array_input(self):
        """Should handle array input."""
        eps0, Mc, a, b = 0.1, 1e12, 0.6, 0.35
        Mh_arr = np.logspace(10, 13, 50)
        
        sfe = star_formation_efficiency(Mh_arr, eps0, Mc, a, b)
        
        assert len(sfe) == len(Mh_arr)
        assert np.all(sfe > 0)
        assert np.all(sfe <= 2 * eps0)


class TestHaloAccretionRate:
    """Tests for halo accretion rate."""
    
    def test_scaling_with_mass(self):
        """Accretion rate should scale with mass."""
        z = 6.0
        M1, M2 = 1e11, 1e12
        
        rate1 = halo_accretion_rate(M1, z)
        rate2 = halo_accretion_rate(M2, z)
        
        # Higher mass should have higher accretion
        assert rate2 > rate1
    
    def test_scaling_with_redshift(self):
        """Accretion rate should increase with redshift."""
        Mh = 1e11
        z1, z2 = 5.0, 7.0
        
        rate1 = halo_accretion_rate(Mh, z1)
        rate2 = halo_accretion_rate(Mh, z2)
        
        # Higher redshift should have higher accretion
        assert rate2 > rate1
    
    def test_positive(self):
        """Accretion rate should always be positive."""
        Mh = np.logspace(9, 13, 20)
        z = 6.0
        
        rates = halo_accretion_rate(Mh, z)
        
        assert np.all(rates > 0)


class TestSFRFromHaloMass:
    """Tests for SFR calculation."""
    
    def test_positive_sfr(self):
        """SFR should be positive."""
        Mh = np.logspace(10, 13, 20)
        z = 6.0
        eps0, Mc, a, b = 0.1, 1e12, 0.6, 0.35
        
        sfr = sfr_from_halo_mass(Mh, z, eps0, Mc, a, b)
        
        assert np.all(sfr > 0)
    
    def test_monotonic_increase(self):
        """SFR should generally increase with mass around Mc."""
        Mh = np.logspace(11, 12, 20)
        z = 6.0
        eps0, Mc, a, b = 0.1, 1e12, 0.6, 0.35
        
        sfr = sfr_from_halo_mass(Mh, z, eps0, Mc, a, b)
        
        # Check monotonicity
        assert np.all(np.diff(sfr) > 0)


class TestUVHaloMassRelation:
    """Tests for UV-halo mass relation."""
    
    def test_initialization(self):
        """Should initialize correctly."""
        uvhmr = UVHaloMassRelation(z=6.0, eps0=0.1, Mc=1e12, a=0.6, b=0.35)
        
        assert uvhmr.z == 6.0
        assert uvhmr.eps0 == 0.1
        assert uvhmr.Mc == 1e12
    
    def test_forward_backward(self):
        """Forward and backward should be inverses."""
        uvhmr = UVHaloMassRelation(z=6.0, eps0=0.1, Mc=1e12, a=0.6, b=0.35)
        
        Mh_orig = 1e11
        MUV = uvhmr.MUV(Mh_orig)
        Mh_recovered = uvhmr.Mhalo(MUV)
        
        # Should recover original mass (within tolerance)
        assert np.isclose(Mh_orig, Mh_recovered, rtol=0.05)
    
    def test_brighter_at_higher_mass(self):
        """Higher mass should give brighter (more negative) MUV."""
        uvhmr = UVHaloMassRelation(z=6.0, eps0=0.1, Mc=1e12, a=0.6, b=0.35)
        
        Mh1, Mh2 = 1e10, 1e12
        MUV1 = uvhmr.MUV(Mh1)
        MUV2 = uvhmr.MUV(Mh2)
        
        # Brighter = more negative
        assert MUV2 < MUV1
    
    def test_sfr_method(self):
        """SFR method should work."""
        uvhmr = UVHaloMassRelation(z=6.0, eps0=0.1, Mc=1e12, a=0.6, b=0.35)
        
        Mh = 1e11
        sfr = uvhmr.sfr(Mh)
        
        assert sfr > 0
        assert isinstance(sfr, (float, np.ndarray))


class TestInverseUVHMR:
    """Tests for inverse UV-halo mass relation."""
    
    def test_inversion(self):
        """Should correctly invert UVHMR."""
        z = 6.0
        eps0, Mc, a, b = 0.1, 1e12, 0.6, 0.35
        
        # Start with a halo mass
        Mh_original = 5e11
        
        # Get MUV
        MUV = MUV_from_halo_mass(Mh_original, z, eps0, Mc, a, b)
        
        # Invert
        Mh_recovered = halo_mass_from_MUV(MUV, z, eps0, Mc, a, b)
        
        # Check recovery
        assert np.isclose(Mh_original, Mh_recovered, rtol=0.05)
    
    def test_array_input(self):
        """Should handle array input."""
        z = 6.0
        eps0, Mc, a, b = 0.1, 1e12, 0.6, 0.35
        
        MUV_arr = np.array([-20, -18, -16])
        Mh_arr = halo_mass_from_MUV(MUV_arr, z, eps0, Mc, a, b)
        
        assert len(Mh_arr) == len(MUV_arr)
        assert np.all(Mh_arr > 0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
