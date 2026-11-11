#!/usr/bin/env python
"""Production fit of the Shuntov+2025b joint clustering + UVLF likelihood.

This is the robust path for a real, long MCMC chain -- run it as a script (NOT
in a notebook), so multiprocessing works cleanly under the ``if __name__ ==
'__main__'`` guard:

    python examples/fit_shuntov25.py --mode fixed_z --bin ha-z5.4 \
        --nwalkers 64 --nsteps 4000 --processes 8 --out chain_z5.4.npz

    python examples/fit_shuntov25.py --mode evolution \
        --nwalkers 64 --nsteps 6000 --processes 16 --out chain_evolution.npz

Each clustering likelihood call is ~0.5-1 s, so parallelise across cores.
``EmceeSampler.run(processes=N)`` builds the (expensive) halomod model once per
worker, so the speed-up is actually realised.

The chain is checkpointed to ``<out>.h5`` every step: re-running the same
command resumes it (delete the file to start afresh). Add ``--converge`` to
stop on the autocorrelation-time criterion instead of after ``--nsteps``.

Caveats of this example relative to the published analysis
----------------------------------------------------------
* The H-alpha/[O III] emitter UVLF points are not in the public data repo (only
  a best-fit model band is). Add them in ``uvlf_points`` to fit the full
  likelihood; until then the fit uses Bouwens+2021 UVLF + clustering.
* The model w(theta) carries no integral-constraint correction, so by default
  only theta < ``--theta-max`` (60 arcsec) is fitted. The two largest-scale
  points of each measurement are consistent with zero or negative, and the
  ~108 arcsec points of the z=5.4 and z=7.3 samples carry quoted errors of
  ~5e-4 that would dominate chi^2 on their own.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "examples"))

from halogal.inference import (  # noqa: E402
    Parameter, ParameterSet, Uniform, Normal,
    ModelParametrization, Dataset, GaussianLikelihood, EmceeSampler,
)
from halogal.config import DEFAULT_REDSHIFT_EVOLUTION  # noqa: E402
from bouwens21_data import bouwens21  # noqa: E402

DATA = REPO / "examples" / "shuntov25_data"

BINS = {
    "ha-z4.3":   dict(z=4.3, bouwens="z4", muv_thresh=-20.0),
    "ha-z5.4":   dict(z=5.4, bouwens="z5", muv_thresh=-20.0),
    "oiii-z7.3": dict(z=7.3, bouwens="z7", muv_thresh=-20.0),
}

# Forwarded to initialize_correlation_model. On the recommended env (numpy<2)
# CAMB works; add transfer_model='EH' if you must run under numpy>=2.
# theta_min/theta_max [arcsec] must bracket the data: the measurements start at
# ~0.4 arcsec, below halogal's default model grid (1 arcsec).
CORR_KW = dict(correlation_type="angular", znum=30, theta_num=40,
               theta_min=0.3, theta_max=3600.0)
HMF_POINTS = 256


def uvlf_points(cfg):
    """UVLF data points: Bouwens+2021 (+ add Ha/[OIII] emitter points here)."""
    bw = bouwens21[cfg["bouwens"]]
    return np.asarray(bw["M_AB"]), np.asarray(bw["Fi_k"]), np.asarray(bw["Fi_k_error"])


def datasets_for(tags, theta_max=60.0):
    out = []
    for tag in tags:
        cfg = BINS[tag]
        z = cfg["z"]
        t = Table.read(DATA / f"w-theta_measurements_{tag}.ecsv")
        theta, w, w_err = (np.array(t[c]) for c in ("theta", "w", "w_err"))
        keep = theta < theta_max
        muv, phi, err = uvlf_points(cfg)
        out.append(Dataset.uvlf(z=z, MUV=muv, phi=phi, sigma=err,
                                label=f"UVLF {tag}"))
        out.append(Dataset.clustering(
            z=z, theta=theta[keep], w=w[keep], sigma=w_err[keep],
            MUV_thresh1=cfg["muv_thresh"],
            corr_kwargs=CORR_KW, label=f"w {tag}"))
    return out


def fixed_z_params():
    return ParameterSet([
        Parameter("eps0", Uniform(0.02, 0.40)),
        Parameter("logMc", Uniform(10.5, 12.5)),
        Parameter("a", Uniform(0.20, 1.20)),
        Parameter("b", Uniform(0.20, 1.50)),
        Parameter("sigma_UV", Uniform(0.20, 1.30)),
        Parameter("logMcut", Uniform(7.0, 12.0)),
        Parameter("logMsat", Uniform(11.0, 14.0)),
        Parameter("asat", Uniform(0.30, 1.50)),
    ])


def evolution_params():
    EV = DEFAULT_REDSHIFT_EVOLUTION
    return ParameterSet([
        Parameter("d_eps0_dz", Normal(EV["d_eps0_dz"], 0.1)),
        Parameter("C_eps0", Uniform(-0.2, 0.6)),
        Parameter("d_logMc_dz", Normal(EV["d_logMc_dz"], 0.2)),
        Parameter("C_logMc", Uniform(9.5, 12.0)),
        Parameter("d_a_dz", Normal(EV["d_a_dz"], 0.1)),
        Parameter("C_a", Uniform(0.2, 1.2)),
        Parameter("d_b_dz", Normal(EV["d_b_dz"], 0.1)),
        Parameter("C_b", Uniform(0.2, 1.5)),
    ])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["fixed_z", "evolution"], default="fixed_z")
    ap.add_argument("--bin", default="ha-z5.4", choices=list(BINS),
                    help="redshift bin for fixed_z mode")
    ap.add_argument("--nwalkers", type=int, default=64)
    ap.add_argument("--nsteps", type=int, default=4000)
    ap.add_argument("--processes", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="chain.npz")
    ap.add_argument("--theta-max", type=float, default=60.0,
                    help="fit w(theta) only below this scale [arcsec]")
    ap.add_argument("--converge", action="store_true",
                    help="stop on the autocorrelation criterion; --nsteps "
                         "becomes the cap")
    args = ap.parse_args()

    if args.mode == "fixed_z":
        params = fixed_z_params()
        data = datasets_for([args.bin], theta_max=args.theta_max)
    else:
        params = evolution_params()
        data = datasets_for(list(BINS), theta_max=args.theta_max)

    like = GaussianLikelihood(
        data, ModelParametrization(args.mode, params), hmf_points=HMF_POINTS)

    # Fail fast if the setup is broken before launching a long chain.
    print("Likelihood sanity check:", like.check())

    print(f"Running emcee: mode={args.mode} nwalkers={args.nwalkers} "
          f"nsteps={args.nsteps} processes={args.processes}")
    checkpoint = str(Path(args.out).with_suffix(".h5"))
    result = EmceeSampler(params, like).run(
        nwalkers=args.nwalkers, nsteps=args.nsteps,
        processes=args.processes, seed=args.seed, progress=True,
        backend=checkpoint, convergence=args.converge)

    print(result)
    print("acceptance fraction:", result.extras["acceptance_fraction"],
          "| autocorrelation time:", result.extras["autocorr_time"])
    best = result.median_params()
    print(f"chi2 at posterior median: {like.chi2(best):.1f} "
          f"for {like.ndata} data points")
    path = result.save(args.out)
    print(f"Saved {len(result.samples)} samples to {path} "
          f"(checkpoint: {checkpoint})")


if __name__ == "__main__":
    main()
