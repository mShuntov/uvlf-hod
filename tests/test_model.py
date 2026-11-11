"""Tests for the forward model: UVHMR, HOD occupations and observables.

Fast unit tests only (no halomod model is built here).
"""

import numpy as np
import pytest

from halogal import HODModel, UVHMRModel
from halogal.config import DEFAULT_HOD_PARAMS
from halogal.cosmology import get_halo_bias, get_halo_mass_function
from halogal.luminosity import MUV_from_SFR, SFR_from_MUV
from halogal.model import Observables


# --------------------------------------------------------------------------
# UV-halo mass relation
# --------------------------------------------------------------------------
class TestUVHMR:
    def test_defaults_come_from_config(self):
        m = UVHMRModel(z=6.0)
        assert m.eps0 == DEFAULT_HOD_PARAMS["eps0"]
        assert m.Mc == pytest.approx(10 ** DEFAULT_HOD_PARAMS["logMc"])

    def test_sfe_double_power_law(self):
        m = UVHMRModel(z=6.0, eps0=0.1, Mc=1e12, a=0.6, b=0.35)
        # 2 eps0 / (x^-a + x^b) equals eps0 at x = 1
        assert m.star_formation_efficiency(1e12)[0] == pytest.approx(0.1)
        low = m.star_formation_efficiency(np.array([1e9, 1e10]))
        assert low[1] / low[0] == pytest.approx(10 ** 0.6, rel=2e-2)
        high = m.star_formation_efficiency(np.array([1e14, 1e15]))
        assert high[1] / high[0] == pytest.approx(10 ** -0.35, rel=2e-2)

    def test_accretion_rate_redshift_scaling(self):
        lo, hi = UVHMRModel(z=4.0), UVHMRModel(z=8.0)
        ratio = hi.halo_accretion_rate(1e11) / lo.halo_accretion_rate(1e11)
        assert ratio[0] == pytest.approx((9.0 / 5.0) ** 2.5)

    def test_sfr_positive_and_monotonic(self):
        m = UVHMRModel(z=6.0)
        sfr = m.sfr(np.logspace(9, 13, 50))
        assert np.all(sfr > 0) and np.all(np.diff(sfr) > 0)

    def test_more_massive_haloes_are_brighter(self):
        m = UVHMRModel(z=6.0)
        muv = m.MUV(np.logspace(9.5, 12.5, 30))
        assert np.all(np.diff(muv) < 0)

    def test_inverse_relation_roundtrip(self):
        m = UVHMRModel(z=6.0)
        Mh = np.array([1e10, 1e11, 1e12])
        # MUV(Mh) is only converged to ~0.02 mag by the dust iteration, which
        # limits the inversion to a few per cent in mass.
        assert m.Mhalo(m.MUV(Mh)) == pytest.approx(Mh, rel=0.1)

    def test_update_parameters(self):
        m = UVHMRModel(z=6.0)
        before = m.MUV(1e11)[0]
        m.update_parameters(eps0=2 * m.eps0)
        assert m.MUV(1e11)[0] < before          # brighter

    def test_sfr_muv_conversion_roundtrip_without_dust(self):
        sfr = np.array([0.1, 1.0, 10.0])
        muv = MUV_from_SFR(sfr, 6.0, add_dust=False)
        assert SFR_from_MUV(muv, 6.0, remove_dust=False) == pytest.approx(sfr)

    def test_dust_only_dims(self):
        sfr = np.array([1.0, 30.0, 300.0])
        dusty = MUV_from_SFR(sfr, 5.0, add_dust=True)
        clean = MUV_from_SFR(sfr, 5.0, add_dust=False)
        assert np.all(dusty >= clean - 1e-9)


# --------------------------------------------------------------------------
# Occupation functions
# --------------------------------------------------------------------------
class TestHOD:
    def test_ncen_is_a_probability_and_monotonic(self):
        m = HODModel(z=6.0)
        Mh = np.logspace(9, 14, 80)
        ncen = m.Ncen(Mh, -19.0)
        assert np.all((ncen >= 0) & (ncen <= 1))
        assert np.all(np.diff(ncen) >= -1e-12)
        assert ncen[0] < 1e-3 and ncen[-1] > 0.99

    def test_brighter_threshold_selects_fewer(self):
        m = HODModel(z=6.0)
        Mh = np.logspace(10, 13, 40)
        assert np.all(m.Ncen(Mh, -21.0) <= m.Ncen(Mh, -19.0))

    def test_nsat_zero_below_mcut(self):
        m = HODModel(z=6.0, Mcut=1e11)
        Mh = np.array([1e9, 5e10, 1e11, 1e12, 1e13])
        nsat = m.Nsat(Mh, -18.0)
        assert np.all(nsat[:3] == 0) and np.all(nsat[3:] > 0)

    def test_integer_masses_are_not_truncated(self):
        m = HODModel(z=6.0)
        as_int = m.Nsat(np.array([10 ** 12, 10 ** 13]), -18.0)
        as_float = m.Nsat(np.array([1e12, 1e13]), -18.0)
        assert as_int == pytest.approx(as_float)

    def test_ngal_is_sum(self):
        m = HODModel(z=6.0)
        Mh = np.logspace(10, 13, 20)
        assert m.Ngal(Mh, -19.0) == pytest.approx(
            m.Ncen(Mh, -19.0) + m.Nsat(Mh, -19.0))


# --------------------------------------------------------------------------
# Halo mass function / bias conventions
# --------------------------------------------------------------------------
class TestCosmology:
    def test_hmf_shape_and_monotonic(self):
        log10m, hmf = get_halo_mass_function(6.0, M_min=8, M_max=13,
                                             num_points=64)
        assert log10m.shape == hmf.shape == (64,)
        assert np.all(hmf > 0) and np.all(np.diff(hmf) < 0)

    def test_hmf_drops_with_redshift_at_high_mass(self):
        _, lo = get_halo_mass_function(5.0, M_min=11, M_max=12, num_points=8)
        _, hi = get_halo_mass_function(8.0, M_min=11, M_max=12, num_points=8)
        assert np.all(hi < lo)

    def test_bias_increases_with_mass(self):
        b = get_halo_bias(np.logspace(10, 13, 10), 6.0)
        assert np.all(np.diff(b) > 0) and b[-1] > 1


# --------------------------------------------------------------------------
# Observables
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def obs():
    o = Observables(HODModel(z=6.0))
    # coarse grid keeps these tests fast; accuracy is checked below
    o._hmf_cache[6.0] = get_halo_mass_function(
        6.0, M_min=6, M_max=15, num_points=256)
    return o


class TestObservables:

    def test_uvlf_positive_and_falls_to_bright_end(self, obs):
        MUV = np.linspace(-22, -17, 11)
        phi = obs.luminosity_function(MUV)
        assert np.all(phi > 0) and np.all(np.diff(phi) > 0)

    def test_integer_magnitudes_give_nonzero_uvlf(self, obs):
        phi_int = obs.luminosity_function(np.arange(-21, -17))
        phi_flt = obs.luminosity_function(np.arange(-21.0, -17.0))
        assert np.all(phi_int > 0) and phi_int == pytest.approx(phi_flt)

    def test_uvlf_integrates_to_number_density(self, obs):
        # n(<M) must equal the integral of Phi over magnitudes brighter than M
        MUV = np.linspace(-26.0, -19.0, 281)
        phi = obs.luminosity_function(MUV)
        n_int = np.trapz(phi, MUV)
        assert n_int == pytest.approx(obs.number_density(-19.0), rel=2e-2)

    def test_inline_parameter_override_persists(self, obs):
        MUV = np.array([-20.0])
        base = obs.luminosity_function(MUV, eps0=0.19)[0]
        boosted = obs.luminosity_function(MUV, eps0=0.30)[0]
        assert boosted > base
        assert obs.hod_model.eps0 == 0.30
        obs.luminosity_function(MUV, eps0=0.19)      # restore for other tests

    def test_unknown_parameter_raises(self, obs):
        with pytest.raises(TypeError, match="unrecognised"):
            obs.luminosity_function(np.array([-20.0]), sigmaUV=0.5)

    def test_bias_grows_with_luminosity(self, obs):
        bg = obs.galaxy_bias(np.array([-21.5, -20.0, -18.5]))
        assert bg[0] > bg[1] > bg[2] > 1

    def test_mean_halo_mass_and_bias_grow_with_luminosity(self, obs):
        assert obs.mean_halo_mass(-21.0) > obs.mean_halo_mass(-19.0)
        assert obs.mean_bias(-21.0) > obs.mean_bias(-19.0)

    def test_coarse_hmf_grid_is_accurate(self):
        MUV = np.linspace(-22, -17, 6)
        fine = Observables(HODModel(z=6.0)).luminosity_function(MUV)
        coarse = Observables(HODModel(z=6.0))
        coarse._hmf_cache[6.0] = get_halo_mass_function(
            6.0, M_min=6, M_max=15, num_points=256)
        assert coarse.luminosity_function(MUV) == pytest.approx(fine, rel=5e-3)
