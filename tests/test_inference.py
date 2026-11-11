"""Tests for the halogal.inference subpackage.

Fast, deterministic checks of the parameter/prior/likelihood machinery, plus
slower end-to-end sampler checks (marked ``slow``). The forward model itself
(colossus/halomod) is exercised only where needed; the clustering interpolation
logic is additionally tested with a lightweight stub so it does not depend on a
working halomod/camb stack.
"""

import numpy as np
import pytest

from halogal.inference import (
    Parameter,
    ParameterSet,
    Uniform,
    Normal,
    LogUniform,
    ModelParametrization,
    Dataset,
    GaussianLikelihood,
    EmceeSampler,
    DynestySampler,
)
from halogal.inference.parametrization import MODEL_PARAM_NAMES


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def _make_uvlf_data(z=5.4, seed=0, frac_err=0.10, n=12):
    """Forward-generate a noisy UVLF from known parameters."""
    from halogal.model import HODModel, Observables

    rng = np.random.default_rng(seed)
    MUV = np.linspace(-22, -17, n)
    obs = Observables(HODModel(z=z))
    truth = dict(eps0=0.19, logMc=11.64)
    mp = ModelParametrization("fixed_z")
    template = Dataset.uvlf(z=z, MUV=MUV, phi=np.ones(n), sigma=np.ones(n))
    phi = template.predict(obs, mp.params_for_redshift(truth, z))
    sigma = frac_err * phi + 1e-6
    phi_obs = phi + sigma * rng.standard_normal(n)
    return MUV, phi_obs, sigma, truth


class _FakeCorrObservables:
    """Stub mimicking Observables.update_correlation_model (no halomod)."""

    def __init__(self, theta_arcsec, w):
        self._sep = np.asarray(theta_arcsec, float)
        self._w = np.asarray(w, float)
        self.received = None

    def update_correlation_model(self, **kwargs):
        self.received = kwargs
        return {"separation": self._sep, "correlation": self._w,
                "type": "angular"}


# --------------------------------------------------------------------------
# Priors
# --------------------------------------------------------------------------
class TestPriors:
    def test_uniform_unit_transform(self):
        u = Uniform(2.0, 6.0)
        assert u.unit_transform(0.0) == pytest.approx(2.0)
        assert u.unit_transform(0.5) == pytest.approx(4.0)
        assert u.unit_transform(1.0) == pytest.approx(6.0)

    def test_uniform_logpdf_support(self):
        u = Uniform(0.0, 2.0)
        assert np.isfinite(u.logpdf(1.0))
        assert u.logpdf(-0.1) == -np.inf
        assert u.logpdf(2.1) == -np.inf

    def test_normal_unit_transform_and_logpdf(self):
        n = Normal(1.0, 2.0)
        assert n.unit_transform(0.5) == pytest.approx(1.0)
        # symmetric quantiles around the mean
        lo = n.unit_transform(0.1586553)
        hi = n.unit_transform(0.8413447)
        assert lo == pytest.approx(1.0 - 2.0, abs=1e-3)
        assert hi == pytest.approx(1.0 + 2.0, abs=1e-3)
        # logpdf maximised at the mean
        assert n.logpdf(1.0) > n.logpdf(3.0)

    def test_loguniform(self):
        lu = LogUniform(1e9, 1e13)
        assert np.log10(lu.unit_transform(0.5)) == pytest.approx(11.0)
        assert lu.logpdf(1e8) == -np.inf
        assert np.isfinite(lu.logpdf(1e11))


# --------------------------------------------------------------------------
# ParameterSet
# --------------------------------------------------------------------------
class TestParameterSet:
    def test_roundtrip_and_fixed(self):
        ps = ParameterSet([
            Parameter("eps0", Uniform(0.0, 1.0)),
            Parameter("logMc", Normal(11.6, 0.5)),
            Parameter("a", value=0.69, fixed=True),
        ])
        assert ps.ndim == 2
        assert ps.names == ["eps0", "logMc"]
        d = ps.to_dict([0.2, 11.7])
        assert d == {"a": 0.69, "eps0": 0.2, "logMc": 11.7}
        # fixed parameter never enters theta
        assert "a" not in ps.names

    def test_log_prior_out_of_bounds(self):
        ps = ParameterSet([Parameter("eps0", Uniform(0.0, 1.0))])
        assert np.isfinite(ps.log_prior([0.5]))
        assert ps.log_prior([1.5]) == -np.inf

    def test_unit_transform_and_initial(self):
        ps = ParameterSet([
            Parameter("eps0", Uniform(0.0, 1.0)),
            Parameter("logMc", Uniform(11.0, 12.0)),
        ])
        theta = ps.unit_transform([0.5, 0.5])
        assert theta == pytest.approx([0.5, 11.5])
        pos = ps.initial(20, rng=np.random.default_rng(0))
        assert pos.shape == (20, 2)
        # every walker must start with finite prior
        assert all(np.isfinite(ps.log_prior(p)) for p in pos)

    def test_fixed_requires_value(self):
        with pytest.raises(ValueError):
            Parameter("a", fixed=True)

    def test_duplicate_names_rejected(self):
        with pytest.raises(ValueError):
            ParameterSet([
                Parameter("eps0", Uniform(0, 1)),
                Parameter("eps0", Uniform(0, 1)),
            ])


# --------------------------------------------------------------------------
# ModelParametrization
# --------------------------------------------------------------------------
class TestParametrization:
    def test_fixed_z_log_to_linear_mass(self):
        ps = ParameterSet([Parameter("logMc", Uniform(11.0, 12.0))])
        mp = ModelParametrization("fixed_z", ps)
        out = mp.params_for_redshift(ps.to_dict([11.7]), z=5.4)
        assert out["Mc"] == pytest.approx(10 ** 11.7)
        # unspecified params filled from config defaults
        assert out["a"] == pytest.approx(0.69)
        assert np.log10(out["Mcut"]) == pytest.approx(9.57)

    def test_evolution_recovers_config_defaults(self):
        # All coefficients fixed at published defaults -> reproduce z=5.4 fit.
        mp = ModelParametrization("evolution")
        out = mp.params_for_redshift({}, z=5.4)
        assert out["eps0"] == pytest.approx(0.198, abs=1e-3)
        assert np.log10(out["Mc"]) == pytest.approx(11.616, abs=1e-3)
        assert out["b"] == pytest.approx(0.646, abs=1e-3)

    def test_unknown_name_rejected(self):
        with pytest.raises(ValueError):
            ModelParametrization(
                "fixed_z", ParameterSet([Parameter("epsilon", Uniform(0, 1))])
            )

    def test_bad_mode(self):
        with pytest.raises(ValueError):
            ModelParametrization("bogus")


# --------------------------------------------------------------------------
# Dataset
# --------------------------------------------------------------------------
class TestDataset:
    def test_sigma_builds_diagonal_cov(self):
        phi = np.array([1.0, 2.0, 3.0])
        sig = np.array([0.1, 0.2, 0.3])
        ds = Dataset.uvlf(z=6, MUV=np.array([-20, -19, -18]), phi=phi, sigma=sig)
        assert np.allclose(ds.cov, np.diag(sig ** 2))

    def test_cholesky_chi2_matches_direct(self):
        phi = np.array([1.0, 2.0, 3.0])
        sig = np.array([0.1, 0.2, 0.3])
        model = np.array([1.1, 1.8, 3.0])
        ds = Dataset.uvlf(z=6, MUV=np.array([-20, -19, -18]), phi=phi, sigma=sig)
        r = phi - model
        assert ds.chi2(model) == pytest.approx(np.sum((r / sig) ** 2))
        # full-matrix covariance gives the identical result
        ds2 = Dataset.uvlf(z=6, MUV=np.array([-20, -19, -18]), phi=phi,
                           cov=np.diag(sig ** 2))
        assert ds2.chi2(model) == pytest.approx(ds.chi2(model))

    def test_requires_exactly_one_of_cov_sigma(self):
        with pytest.raises(ValueError):
            Dataset.uvlf(z=6, MUV=[-20, -19], phi=[1, 2])
        with pytest.raises(ValueError):
            Dataset.uvlf(z=6, MUV=[-20, -19], phi=[1, 2],
                         sigma=[0.1, 0.1], cov=np.eye(2))

    def test_clustering_requires_threshold(self):
        with pytest.raises(ValueError):
            Dataset.clustering(z=6, theta=[1, 2], w=[1, 2], sigma=[1, 1])

    def test_clustering_predict_interpolation_with_stub(self):
        # log-theta interpolation of model w onto the data abscissa.
        theta_model = np.logspace(0, 3.5, 50)       # 1 .. ~3000 arcsec
        w_model = 5.0 * (theta_model / 10.0) ** -0.8
        fake = _FakeCorrObservables(theta_model, w_model)
        theta_data = np.array([2.0, 20.0, 200.0])
        ds = Dataset.clustering(z=6, theta=theta_data, w=np.ones(3),
                                sigma=np.ones(3), MUV_thresh1=-19.0)
        model_params = {n: 1.0 for n in MODEL_PARAM_NAMES}
        model_params["junk"] = 999.0  # must be filtered out
        pred = ds.predict(fake, model_params)
        expected = np.interp(np.log10(theta_data), np.log10(theta_model), w_model)
        assert pred == pytest.approx(expected)
        # only recognised correlation params are forwarded
        assert "junk" not in fake.received
        assert set(fake.received) <= set(MODEL_PARAM_NAMES)


# --------------------------------------------------------------------------
# GaussianLikelihood (UVLF only; no halomod needed)
# --------------------------------------------------------------------------
class TestGaussianLikelihoodUVLF:
    def test_likelihood_peaks_near_truth(self):
        MUV, phi, sigma, truth = _make_uvlf_data()
        params = ParameterSet([
            Parameter("eps0", Uniform(0.05, 0.4)),
            Parameter("logMc", Uniform(11.0, 12.3)),
        ])
        like = GaussianLikelihood(
            [Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=sigma)],
            ModelParametrization("fixed_z", params), catch_errors=False)
        ll_truth = like.log_likelihood({"eps0": 0.19, "logMc": 11.64})
        ll_bad = like.log_likelihood({"eps0": 0.35, "logMc": 11.1})
        assert ll_truth > ll_bad

    def test_too_many_clustering_per_bin_rejected(self):
        d1 = Dataset.clustering(z=6, theta=[1, 2], w=[1, 1], sigma=[1, 1],
                                MUV_thresh1=-19, label="w bright")
        d2 = Dataset.clustering(z=6, theta=[1, 2], w=[1, 1], sigma=[1, 1],
                                MUV_thresh1=-20, label="w faint")
        with pytest.raises(ValueError, match="one clustering dataset"):
            GaussianLikelihood([d1, d2], ModelParametrization("fixed_z"))

    def test_hmf_points_speeds_uvlf_without_changing_result(self):
        # A coarser HMF grid for a clustering-free bin must barely change the
        # (smooth) UVLF likelihood.
        MUV, phi, sigma, _ = _make_uvlf_data()
        ds = Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=sigma)
        full = GaussianLikelihood([ds], ModelParametrization("fixed_z"))
        coarse = GaussianLikelihood([ds], ModelParametrization("fixed_z"),
                                    hmf_points=256)
        ll_full = full.log_likelihood({"eps0": 0.19})
        ll_coarse = coarse.log_likelihood({"eps0": 0.19})
        assert abs(ll_coarse - ll_full) < 0.01 * abs(ll_full)
        # the coarse bin actually used the 256-point grid
        coarse._ensure_setup()
        log10m, _ = coarse._bins[5.4]["obs"]._get_hmf()
        assert len(log10m) == 256

    def test_check_returns_finite_per_dataset(self):
        MUV, phi, sigma, _ = _make_uvlf_data()
        like = GaussianLikelihood(
            [Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=sigma, label="uvlf")],
            ModelParametrization("fixed_z"))
        out = like.check({"eps0": 0.19})
        assert set(out) == {"uvlf"} and np.isfinite(out["uvlf"])

    def test_setup_error_fails_fast_even_with_catch_errors(self):
        # A clustering dataset with a bogus halomod model makes
        # initialize_correlation_model raise; that is a configuration error and
        # must propagate from log_likelihood even though catch_errors=True
        # (it must NOT be silently masked as -inf).
        pytest.importorskip("halomod")
        theta = np.array([1.0, 10.0, 100.0])
        clu = Dataset.clustering(
            z=5.4, theta=theta, w=np.ones(3), sigma=np.ones(3), MUV_thresh1=-20.0,
            corr_kwargs=dict(correlation_type="angular",
                             exclusion_model="NOT_A_REAL_MODEL"))
        like = GaussianLikelihood([clu], ModelParametrization("fixed_z"),
                                  catch_errors=True)
        with pytest.raises(Exception):
            like.log_likelihood({"eps0": 0.19})

    def test_pickle_drops_live_observables(self):
        import pickle
        MUV, phi, sigma, _ = _make_uvlf_data()
        like = GaussianLikelihood(
            [Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=sigma)],
            ModelParametrization("fixed_z"))
        like.log_likelihood({"eps0": 0.19})        # build live bins
        restored = pickle.loads(pickle.dumps(like))
        assert restored._bins is None              # rebuilt lazily per process
        assert np.isfinite(restored.log_likelihood({"eps0": 0.19}))


# --------------------------------------------------------------------------
# Samplers (slow)
# --------------------------------------------------------------------------
class TestSamplers:
    @pytest.mark.slow
    def test_emcee_uvlf_recovery(self):
        MUV, phi, sigma, truth = _make_uvlf_data(seed=1)
        params = ParameterSet([
            Parameter("eps0", Uniform(0.05, 0.4)),
            Parameter("logMc", Uniform(11.0, 12.3)),
        ])
        like = GaussianLikelihood(
            [Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=sigma)],
            ModelParametrization("fixed_z", params), catch_errors=False,
            hmf_points=256)
        res = EmceeSampler(params, like).run(
            nwalkers=20, nsteps=300, progress=False, seed=1)
        s = res.summary()
        assert abs(s["eps0"]["median"] - 0.19) < 0.06
        assert abs(s["logMc"]["median"] - 11.64) < 0.25

    def test_processes_and_pool_are_mutually_exclusive(self):
        MUV, phi, sigma, _ = _make_uvlf_data()
        params = ParameterSet([Parameter("eps0", Uniform(0.05, 0.4))])
        like = GaussianLikelihood(
            [Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=sigma)],
            ModelParametrization("fixed_z", params))
        with pytest.raises(ValueError):
            EmceeSampler(params, like).run(processes=2, pool=object())

    @pytest.mark.slow
    def test_emcee_processes_parallel_runs(self):
        # The processes= path builds the likelihood once per worker and returns
        # finite results (UVLF-only keeps the per-call cost low).
        MUV, phi, sigma, _ = _make_uvlf_data(seed=2)
        params = ParameterSet([
            Parameter("eps0", Uniform(0.05, 0.4)),
            Parameter("logMc", Uniform(11.0, 12.3)),
        ])
        like = GaussianLikelihood(
            [Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=sigma)],
            ModelParametrization("fixed_z", params))
        res = EmceeSampler(params, like).run(
            nwalkers=8, nsteps=10, processes=2, progress=False, seed=0)
        assert res.samples.shape[1] == 2
        assert np.all(np.isfinite(res.extras["log_prob"]))

    @pytest.mark.slow
    def test_dynesty_uvlf_recovery_and_evidence(self):
        MUV, phi, sigma, truth = _make_uvlf_data(seed=0)
        params = ParameterSet([
            Parameter("eps0", Uniform(0.05, 0.4)),
            Parameter("logMc", Uniform(11.0, 12.3)),
        ])
        like = GaussianLikelihood(
            [Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=sigma)],
            ModelParametrization("fixed_z", params), catch_errors=False,
            hmf_points=256)
        res = DynestySampler(params, like).run(
            nlive=100, dlogz=0.1, print_progress=False, seed=2)
        s = res.summary()
        assert np.isfinite(res.log_evidence)
        assert abs(s["eps0"]["median"] - 0.19) < 0.06
        assert abs(s["logMc"]["median"] - 11.64) < 0.25


# --------------------------------------------------------------------------
# Live clustering forward model (slow; skips on env incompatibility)
# --------------------------------------------------------------------------
class TestClusteringLive:
    @pytest.mark.slow
    def test_clustering_predict_finite(self):
        pytest.importorskip("halomod")
        theta = np.array([1.5, 5.0, 20.0, 100.0, 500.0, 2000.0])
        clu = Dataset.clustering(
            z=5.4, theta=theta, w=np.ones_like(theta),
            sigma=0.2 * np.ones_like(theta), MUV_thresh1=-18.0,
            corr_kwargs=dict(znum=30, theta_num=25))
        like = GaussianLikelihood([clu], ModelParametrization("fixed_z"),
                                  catch_errors=False)
        try:
            pred = like.predict_all({"eps0": 0.19})[clu.label]["model"]
        except Exception as exc:  # halomod/camb/version incompatibility
            pytest.skip(f"halomod model construction failed in this env: {exc}")
        assert pred.shape == theta.shape
        assert np.all(np.isfinite(pred))


# --------------------------------------------------------------------------
# Convergence loop, HDF backend/resume, and save/load
# --------------------------------------------------------------------------
class TestConvergenceAndCheckpoint:
    def _tiny_problem(self, seed=0):
        MUV, phi, sigma, truth = _make_uvlf_data(seed=seed, n=8)
        params = ParameterSet([
            Parameter("eps0", Uniform(0.05, 0.4)),
            Parameter("logMc", Uniform(11.0, 12.3)),
        ])
        like = GaussianLikelihood(
            [Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=sigma)],
            ModelParametrization("fixed_z", params),
            catch_errors=False, hmf_points=256)
        return params, like, truth

    @pytest.mark.slow
    def test_convergence_stops_early(self):
        params, like, _ = self._tiny_problem(seed=1)
        res = EmceeSampler(params, like).run(
            nwalkers=16, convergence=True, max_steps=4000, check_every=100,
            n_tau=25, tau_change=0.05, progress=False, seed=1)
        assert res.extras["converged"] is True
        # stopped before the cap
        assert res.extras["n_iterations"] < 4000
        assert res.extras["autocorr_time"] is not None

    @pytest.mark.slow
    def test_max_steps_cap_when_not_converged(self):
        params, like, _ = self._tiny_problem(seed=2)
        # absurd criterion (n_tau huge) -> never "converged" -> hits the cap
        res = EmceeSampler(params, like).run(
            nwalkers=16, convergence=True, max_steps=300, check_every=100,
            n_tau=100000, progress=False, seed=2)
        assert res.extras["converged"] is False
        assert res.extras["n_iterations"] == 300

    @pytest.mark.slow
    def test_hdf_backend_resume_grows_chain(self, tmp_path):
        params, like, _ = self._tiny_problem(seed=3)
        h5 = str(tmp_path / "chain.h5")
        EmceeSampler(params, like).run(
            nwalkers=16, nsteps=100, backend=h5, progress=False, seed=3)
        import emcee
        assert emcee.backends.HDFBackend(h5).iteration == 100
        # resume: same backend, more steps, iteration count grows
        EmceeSampler(params, like).run(
            nwalkers=16, nsteps=50, backend=h5, progress=False, seed=4)
        assert emcee.backends.HDFBackend(h5).iteration == 150

    @pytest.mark.slow
    def test_convergence_resume_from_backend(self, tmp_path):
        params, like, _ = self._tiny_problem(seed=5)
        h5 = str(tmp_path / "chain.h5")
        EmceeSampler(params, like).run(
            nwalkers=16, nsteps=100, backend=h5, progress=False, seed=5)
        # convergence-mode resume must read the last state from the backend
        res = EmceeSampler(params, like).run(
            nwalkers=16, convergence=True, max_steps=100, check_every=50,
            n_tau=100000, backend=h5, progress=False, seed=6)
        assert res.extras["n_iterations"] == 200

    def test_ndim_mismatch_on_resume_raises(self, tmp_path):
        params, like, _ = self._tiny_problem(seed=7)
        h5 = str(tmp_path / "chain.h5")
        EmceeSampler(params, like).run(
            nwalkers=16, nsteps=20, backend=h5, progress=False, seed=7)
        # a different-dimension parameter set cannot resume the same backend
        p2 = ParameterSet([Parameter("eps0", Uniform(0.05, 0.4))])
        like2 = GaussianLikelihood(
            [Dataset.uvlf(z=5.4, MUV=np.linspace(-22, -17, 8),
                          phi=np.ones(8), sigma=np.ones(8))],
            ModelParametrization("fixed_z", p2), hmf_points=256)
        with pytest.raises(ValueError):
            EmceeSampler(p2, like2).run(
                nwalkers=16, nsteps=20, backend=h5, progress=False)

    def test_result_save_load_roundtrip(self, tmp_path):
        from halogal.inference.results import InferenceResult
        names = ["a", "b"]
        samples = np.random.default_rng(0).normal(size=(200, 2))
        res = InferenceResult(
            names=names, samples=samples, backend="emcee",
            extras={"converged": True, "burnin": 33, "thin": 1,
                    "acceptance_fraction": 0.4,
                    "autocorr_time": np.array([12.0, 15.0]),
                    "sampler": object()})   # live object must be dropped
        path = res.save(str(tmp_path / "out"))
        loaded = InferenceResult.load(path)
        assert loaded.names == names
        assert np.allclose(loaded.samples, samples)
        assert loaded.extras["converged"] is True
        assert loaded.extras["burnin"] == 33
        assert "sampler" not in loaded.extras
        np.testing.assert_allclose(
            loaded.extras["autocorr_time"], [12.0, 15.0])


# --------------------------------------------------------------------------
# Input validation and robustness guards
# --------------------------------------------------------------------------
class TestValidation:
    def test_dataset_rejects_bad_inputs(self):
        with pytest.raises(ValueError, match="non-finite"):
            Dataset.uvlf(z=6, MUV=[-20, -19], phi=[1e-3, np.nan],
                         sigma=[1e-4, 1e-4])
        with pytest.raises(ValueError, match="strictly positive"):
            Dataset.uvlf(z=6, MUV=[-20, -19], phi=[1e-3, 1e-3],
                         sigma=[1e-4, 0.0])
        with pytest.raises(ValueError, match="symmetric"):
            Dataset.uvlf(z=6, MUV=[-20, -19], phi=[1e-3, 1e-3],
                         cov=[[1.0, 0.5], [0.0, 1.0]])
        with pytest.raises(ValueError, match="positive definite"):
            Dataset.uvlf(z=6, MUV=[-20, -19], phi=[1e-3, 1e-3],
                         cov=[[1.0, 2.0], [2.0, 1.0]])

    def test_clustering_threshold_order_and_theta(self):
        with pytest.raises(ValueError, match="brighter"):
            Dataset.clustering(z=6, theta=[1, 2], w=[1, 1], sigma=[1, 1],
                               MUV_thresh1=-20, MUV_thresh2=-19)
        with pytest.raises(ValueError, match="strictly positive"):
            Dataset.clustering(z=6, theta=[0, 2], w=[1, 1], sigma=[1, 1],
                               MUV_thresh1=-20)

    def test_separation_range_is_enforced(self):
        clu = Dataset.clustering(z=6, theta=[0.5, 10.0], w=[1, 1],
                                 sigma=[1, 1], MUV_thresh1=-20)
        clu.check_separation_range(np.logspace(-1, 3, 20))     # inside: ok
        with pytest.raises(ValueError, match="model grid"):
            clu.check_separation_range(np.logspace(0, 3, 20))  # 0.5 < 1.0

    def test_duplicate_labels_rejected(self):
        MUV, phi, sigma, _ = _make_uvlf_data()
        d1 = Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=sigma)
        d2 = Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=sigma)
        with pytest.raises(ValueError, match="unique"):
            GaussianLikelihood([d1, d2], ModelParametrization("fixed_z"))

    def test_fixed_z_with_several_redshifts_warns(self):
        MUV, phi, sigma, _ = _make_uvlf_data()
        d1 = Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=sigma)
        d2 = Dataset.uvlf(z=6.0, MUV=MUV, phi=phi, sigma=sigma)
        with pytest.warns(UserWarning, match="fixed_z"):
            GaussianLikelihood([d1, d2], ModelParametrization("fixed_z"))

    def test_initial_value_outside_prior_rejected(self):
        with pytest.raises(ValueError, match="outside its prior"):
            Parameter("eps0", Uniform(0.05, 0.4), value=0.9)

    def test_unphysical_evolved_parameters_give_minus_inf(self):
        # sigma_UV(z) = -0.2 z + 0.74 < 0 at z = 5.4. Without the guard the
        # model returns a finite (meaningless) UVLF and a finite likelihood.
        MUV, phi, sigma, _ = _make_uvlf_data()
        like = GaussianLikelihood(
            [Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=sigma, label="u")],
            ModelParametrization("evolution"), catch_errors=False,
            hmf_points=256)
        assert np.isfinite(like.log_likelihood({}))
        bad = {"d_sigmaUV_dz": -0.2}
        assert like.log_likelihood(bad) == -np.inf
        assert like.chi2(bad) == np.inf
        with pytest.raises(ValueError, match="Unphysical"):
            like.check(bad)
        assert like.ndata == MUV.size

    def test_is_physical(self):
        mp = ModelParametrization("fixed_z")
        good = mp.params_for_redshift({}, 5.4)
        assert mp.is_physical(good)
        assert not mp.is_physical(dict(good, eps0=-0.1))
        assert not mp.is_physical(dict(good, sigma_UV=0.0))
        assert not mp.is_physical(dict(good, a=np.nan))
        # an absurd log-mass overflows to inf instead of raising
        assert not mp.is_physical(mp.params_for_redshift({"logMc": 400.0}, 5.4))

    def test_bad_initial_walkers_rejected(self):
        MUV, phi, sigma, _ = _make_uvlf_data()
        params = ParameterSet([Parameter("eps0", Uniform(0.05, 0.4))])
        like = GaussianLikelihood(
            [Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=sigma)],
            ModelParametrization("fixed_z", params), hmf_points=256)
        sampler = EmceeSampler(params, like)
        with pytest.raises(ValueError, match="shape"):
            sampler.run(nwalkers=8, nsteps=2, initial=np.full((4, 1), 0.2),
                        progress=False)
        with pytest.raises(ValueError, match="outside the prior"):
            sampler.run(nwalkers=8, nsteps=2, initial=np.full((8, 1), 0.9),
                        progress=False)

    def test_setup_error_surfaces_before_workers_start(self):
        # With processes=, a setup error raised only inside the Pool
        # initializer makes multiprocessing respawn workers forever (a silent
        # hang). It must be raised in the parent instead.
        pytest.importorskip("halomod")
        theta = np.array([2.0, 10.0, 100.0])
        clu = Dataset.clustering(
            z=5.4, theta=theta, w=np.ones(3), sigma=np.ones(3),
            MUV_thresh1=-20.0,
            corr_kwargs=dict(exclusion_model="NOT_A_REAL_MODEL"))
        params = ParameterSet([Parameter("eps0", Uniform(0.05, 0.4))])
        like = GaussianLikelihood([clu], ModelParametrization("fixed_z", params))
        with pytest.raises(Exception):
            EmceeSampler(params, like).run(
                nwalkers=8, nsteps=2, processes=2, progress=False)
        with pytest.raises(Exception):
            DynestySampler(params, like).run(
                nlive=20, processes=2, print_progress=False)

    def test_dynesty_processes_and_pool_are_mutually_exclusive(self):
        MUV, phi, sigma, _ = _make_uvlf_data()
        params = ParameterSet([Parameter("eps0", Uniform(0.05, 0.4))])
        like = GaussianLikelihood(
            [Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=sigma)],
            ModelParametrization("fixed_z", params), hmf_points=256)
        with pytest.raises(ValueError):
            DynestySampler(params, like).run(processes=2, pool=object())


class TestSamplerAgreement:
    @pytest.mark.slow
    def test_dynesty_processes_matches_emcee(self):
        # Same problem through both backends (dynesty in parallel): the
        # posteriors must agree to well within their widths.
        MUV, phi, sigma, _ = _make_uvlf_data(seed=4)
        params = ParameterSet([
            Parameter("eps0", Uniform(0.05, 0.4)),
            Parameter("logMc", Uniform(11.0, 12.3)),
        ])
        like = GaussianLikelihood(
            [Dataset.uvlf(z=5.4, MUV=MUV, phi=phi, sigma=sigma)],
            ModelParametrization("fixed_z", params), catch_errors=False,
            hmf_points=256)
        em = EmceeSampler(params, like).run(
            nwalkers=20, nsteps=250, progress=False, seed=4).summary()
        dy = DynestySampler(params, like).run(
            nlive=100, dlogz=0.1, processes=2, print_progress=False,
            seed=4).summary()
        for name in ("eps0", "logMc"):
            width = 0.5 * (em[name]["upper"] + em[name]["lower"])
            assert abs(em[name]["median"] - dy[name]["median"]) < width
            ratio = (dy[name]["upper"] + dy[name]["lower"]) / (2 * width)
            assert 0.5 < ratio < 2.0
