"""MCMC inference for the ``halogal`` forward model.

This subpackage fits the UVHMR/HOD model to data with gradient-free samplers
(emcee and dynesty), behind one interface. The user controls **what data** to
fit (any mix of UVLF, clustering, and number-density datasets), **which
parameters** are free (per-parameter ``fixed`` flag), **which redshift regime**
(single fixed-z bin or joint redshift evolution across bins), and **which
sampler**.

Layers
------
``parameters``      Priors and the ``ParameterSet`` (theta <-> named-dict map).
``parametrization`` ``ModelParametrization``: theta -> per-redshift model params.
``data``            ``Dataset`` containers (covariance-aware).
``likelihood``      ``GaussianLikelihood`` composed from a list of datasets.
``samplers``        ``EmceeSampler`` and ``DynestySampler``.
``results``         ``InferenceResult`` (summary, evidence, corner plot).

Example
-------
>>> from halogal.inference import (
...     Dataset, GaussianLikelihood, ParameterSet, Parameter,
...     Uniform, Normal, ModelParametrization, EmceeSampler)
>>> data = [Dataset.uvlf(z=6, MUV=MUV, phi=phi, sigma=phi_err),
...         Dataset.clustering(z=6, theta=th, w=w, cov=w_cov,
...                            MUV_thresh1=-19.1)]
>>> params = ParameterSet([
...     Parameter("eps0", Uniform(0.01, 0.5)),
...     Parameter("logMc", Normal(11.6, 0.5)),
...     Parameter("a", Uniform(0.3, 1.2)),
...     Parameter("b", value=0.65, fixed=True),
...     Parameter("sigma_UV", Uniform(0.1, 1.5)),
...     Parameter("logMcut", Uniform(7.0, 12.0)),
...     Parameter("logMsat", Uniform(11.0, 14.0)),
...     Parameter("asat", Uniform(0.2, 1.5))])
>>> like = GaussianLikelihood(
...     data, ModelParametrization("fixed_z", params))
>>> result = EmceeSampler(params, like).run(nwalkers=64, nsteps=5000)
>>> print(result.summary())
"""

from .parameters import (
    Prior,
    Uniform,
    Normal,
    LogUniform,
    Parameter,
    ParameterSet,
)
from .parametrization import (
    ModelParametrization,
    HOD_PARAM_SPEC,
    MODEL_PARAM_NAMES,
    fixed_z_param_names,
    evolution_param_names,
)
from .data import Dataset
from .likelihood import GaussianLikelihood
from .samplers import EmceeSampler, DynestySampler
from .results import InferenceResult

__all__ = [
    "Prior",
    "Uniform",
    "Normal",
    "LogUniform",
    "Parameter",
    "ParameterSet",
    "ModelParametrization",
    "HOD_PARAM_SPEC",
    "MODEL_PARAM_NAMES",
    "fixed_z_param_names",
    "evolution_param_names",
    "Dataset",
    "GaussianLikelihood",
    "EmceeSampler",
    "DynestySampler",
    "InferenceResult",
]
