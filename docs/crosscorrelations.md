# Cross-correlations

`halogal.crosscorr` computes the halo-model cross-correlation between two
UV-magnitude-selected galaxy samples drawn from the same underlying HOD.

The motivating case is a sparse, bright, spectroscopically confirmed sample
(e.g. Euclid $z>7$ dropouts confirmed with JWST/NIRCam WFSS) cross-correlated
against a fainter, denser sample found in the same fields (e.g. [O III]
emitters in the grism, or a photometric catalogue). The cross-correlation
measures which haloes host the bright galaxies — information the UV luminosity
function alone is degenerate against.

## Why cross-correlate

At fixed halo mass the HOD places a galaxy's UV magnitude with scatter
$\sigma_{\rm UV}$. Increasing $\sigma_{\rm UV}$ up-scatters abundant low-mass
haloes into bright magnitude bins, which can be compensated by lowering the
star-formation efficiency $\epsilon_0$ and shifting the characteristic mass
$M_{\rm c}$. The resulting models can reproduce the same luminosity function
to far better than the measurement errors while populating bright bins with
*completely different haloes*. Clustering — and in particular the number of
faint companions sharing a halo with a bright galaxy — breaks that degeneracy.

## Sample definitions

A `Selection` picks galaxies with

$$M_{\rm UV}^{\rm bright} < M_{\rm UV} < M_{\rm UV}^{\rm faint}.$$

UV absolute magnitudes are negative and brighter objects are more negative, so
`bright` is the more negative bound. Omitting it gives an open-ended threshold
sample.

Occupations for a bin follow by differencing two threshold occupations:

$$N_c^X(M) = N_c(M; M_{\rm UV}^{\rm faint}) - N_c(M; M_{\rm UV}^{\rm bright}),$$

with halogal's error-function central occupation

$$N_c(M; M_{\rm UV}^{\rm thr}) = \frac{1}{2}\left[1 - \mathrm{erf}\!\left(
\frac{\mathcal{M}(M) - M_{\rm UV}^{\rm thr}}{\sqrt{2}\,\sigma_{\rm UV}}\right)\right],$$

where $\mathcal{M}(M)$ is the mean UV magnitude at halo mass $M$, and the
satellite occupation

$$N_s^X(M) = N_c^X(M)\left(\frac{M - M_{\rm cut}}{M_{\rm sat}}\right)^{\alpha_{\rm sat}}
\Theta(M - M_{\rm cut}).$$

The mean comoving number density is

$$\bar n_X = \int d\ln M\, \frac{dn}{d\ln M}\left[N_c^X(M) + N_s^X(M)\right].$$

**The two selections should be disjoint.** If they overlap, a galaxy belongs to
both samples and contributes a self-pair at zero separation, which is not part
of the correlation function. `CrossCorrelation` warns when it detects overlap.

## One-halo term

$$P^{1h}_{12}(k) = \frac{1}{\bar n_1 \bar n_2}\int d\ln M\, \frac{dn}{d\ln M}
\Big[\langle N_c^1 N_s^2\rangle\, u(k|M)
   + \langle N_s^1 N_c^2\rangle\, u(k|M)
   + \langle N_s^1 N_s^2\rangle\, u^2(k|M)\Big]$$

where $u(k|M)$ is the normalised Fourier transform of the halo profile (NFW by
default), taken from the halomod backbone.

**There is no central–central term.** A halo hosts at most one central, so for
disjoint selections that central belongs to at most one sample; a
central–central "pair" would be the same object counted twice.

The physical content of the measurement is **the mean number of faint
companions sharing a halo with a bright galaxy.** $\sigma_{\rm UV}$ controls
which halo mass hosts the bright galaxy, and therefore how many satellites
accompany it.

### Do not drop the host-satellite terms

It is tempting to argue that the bright sample is rare and therefore
essentially all centrals ($N_s^1 \approx 0$), leaving only the $cs$ term.
**That approximation is wrong by a factor of two, and `CrossCorrelation` does
not make it** — `power_1h` always evaluates all three terms.

The reason is structural. `halogal` parametrises satellites as

$$N_s^X(M) = N_c^X(M)\, f(M), \qquad
f(M) = \left[\frac{M - M_{\rm cut}}{M_{\rm sat}}\right]^{\alpha_{\rm sat}},$$

with the *same* $f(M)$ for every selection. Hence

$$\langle N_c^1 N_s^2\rangle = N_c^1 N_c^2 f, \qquad
\langle N_s^1 N_c^2\rangle = N_c^1 f\, N_c^2,$$

which are **identically equal**, however rare sample 1 is. Every host-emitter
pair requires exactly one satellite; there are two equally likely choices of
which member it is. Dropping $N_s^1$ discards one of them.

Measured at $z = 8$ for hosts $M_{\rm UV} < -22$ crossed with emitters
$-22 < M_{\rm UV} < -20$ (fiducial parameters), the one-halo pair budget is

| term | contribution |
|---|---|
| $cs$ (host central $\times$ emitter satellite) | 48.4% |
| $sc$ (host satellite $\times$ emitter central) | 48.4% |
| $ss$ (host satellite $\times$ emitter satellite) | 3.3% |

so setting $N_s^1 = 0$ loses 51.6% of the one-halo pairs. In $w_p$ that is a
factor 1.87 at $r_p \lesssim 0.2$ cMpc, falling to 1.002 beyond
$r_p \simeq 0.35$ cMpc where the two-halo term takes over. The two-halo leg is
insensitive to the same question: host satellites carry $b = 12.70$ against
$b = 12.07$ for host centrals, so removing them shifts large-scale $w_p$ by
0.36%.

**A related consequence of this parametrization.** Because $N_s = N_c f(M)$
carries no requirement that the central outshine its satellites, the satellite
fraction *rises* toward brighter selections — 2.6% at $M_{\rm UV} < -20$, 7.0%
at $-22$, 11.1% at $-23$ — since brighter samples occupy more massive haloes
where $f(M)$ is larger. This is the opposite of the usual expectation that the
brightest galaxies are overwhelmingly centrals. A model placing satellites in
subhaloes (lower mass, hence fainter) would suppress bright satellites. Treat
the bright-end satellite fraction as an assumption of the occupation form, not
as a derived result.

### Correlated scatter

The pair counts allow the two samples' occupations to be correlated:

$$\langle N_c^1 N_s^2\rangle = N_c^1 N_s^2 + R_{cs}\,\sigma_c^1 \sigma_s^2, \qquad
\langle N_s^1 N_s^2\rangle = N_s^1 N_s^2 + R_{ss}\,\sigma_s^1 \sigma_s^2,$$

with Bernoulli centrals $\sigma_c = \sqrt{N_c(1 - N_c)}$ and Poisson satellites
$\sigma_s = \sqrt{N_s}$.

The default $R_{cs} = R_{sc} = R_{ss} = 0$ assumes the $M_{\rm UV}$ scatter is
**independent between a central and its satellites**. This is an assumption,
not a derivation. If scatter is partly environmental or assembly-driven then it
is correlated within a halo and $R_{cs} \neq 0$ — the same physics as assembly
bias. Treat $R_{cs}$ as a nuisance parameter rather than asserting zero.

**The attainable range of $R$ is much smaller than order unity.** A pair count
cannot be negative, so $R_{cs} \gtrsim -N_c^1 N_s^2/(\sigma_c^1\sigma_s^2)$.
For a rare bright selection $N_c \ll 1$ while $\sigma_c=\sqrt{N_c(1-N_c)}
\approx\sqrt{N_c}$, so the bound is roughly $|R| \lesssim \sqrt{N_c N_s}$,
which can be $10^{-2}$ or smaller. Values of order 0.3 drive the one-halo term
negative; `CrossCorrelation` warns and clips when that happens. Scan $R$ over a
range justified by the selection's own occupations, not a fixed 0.3.

## Two-halo term

$$P^{2h}_{12}(k) = b_1(k)\, b_2(k)\, P_{hh}(k), \qquad
b_X(k) = \frac{1}{\bar n_X}\int d\ln M\, \frac{dn}{d\ln M}\, b(M)
\left[N_c^X + N_s^X u(k|M)\right]$$

where $P_{hh}$ is the halo-centre power spectrum from the backbone (nonlinear
by default). As $k \to 0$, $u \to 1$ and $b_X(k) \to b_X$, the number-weighted
effective bias returned by `bias_effective`, so

$$P^{2h}_{12}(k) \to b_1 b_2 P_{hh}(k).$$

This limit is a useful correctness check and is verified in the test suite.

## Real space and projections

$$\xi_{12}(r) = \frac{1}{2\pi^2}\int dk\, k^2\, P_{12}(k)\, j_0(kr),
\qquad P_{12} = P^{1h}_{12} + P^{2h}_{12}$$

evaluated by Hankel transform (Ogata's method, via halomod's utilities).

### Projected, $w_p(r_p)$

$$w_p(r_p) = 2\int_0^{\pi_{\rm max}} d\pi\,
\xi_{12}\!\left(\sqrt{r_p^2 + \pi^2}\right)$$

Appropriate when both samples have spectroscopic redshifts, so line-of-sight
separations are known. This is the preferred observable when available: it
retains the line-of-sight information that projection would otherwise dilute.

### Angular, $w(\theta)$

$$w_{12}(\theta) = \int d\chi\, p_1(\chi)\, p_2(\chi)
\int du\, \xi_{12}\!\left(\sqrt{u^2 + \chi^2\theta^2}\right)$$

in the Limber approximation, with the redshift distributions normalised such
that $\int p_i(\chi)\, d\chi = 1$. Needed when one sample is photometric.

Cross-correlating a narrow spectroscopic distribution against a broad
photometric one dilutes the signal roughly by the fraction of the photometric
sample lying inside the spectroscopic window — captured by the
$\int p_1 p_2\, d\chi$ overlap integral. This is the trade-off between a small,
clean spectroscopic sample and a large, diluted photometric one.

The correlation function is evaluated at the model's own redshift rather than
evolved across the window. That is fine for a narrow selection and approximate
for a broad one.

### Stacked counts

With only a handful of bright hosts, a correlation-function estimator is
Poisson-starved and an integral constraint bites. Counting companions in
cylinders is more robust:

$$\langle N \rangle = \bar n_2 \left[\pi r_p^2 \cdot 2\pi_{\rm max}
+ 2\pi \int_0^{r_p} R\, w_p(R)\, dR\right]$$

`mean_companions` returns this alongside the unclustered expectation, so the
forecast reduces to Poisson statistics on an integer count.

## Implementation notes

All ingredients — halo mass function, halo bias, tracer profile, halo-centre
power — come from **a single halomod `TracerHaloModel` backbone**, shared by
both samples. This is deliberate. halogal has historically been able to serve
the HMF from two different code paths whose results differ by up to 0.47 dex;
using one backbone makes that class of error impossible here.

`dlog10m` defaults to **0.02**, not halomod's 0.1. The central occupation is an
error function of width $\sigma_{\rm UV}$ in magnitude; translated to halo mass
this can be narrower than one grid cell of a coarse mass grid, in which case
the mass integrals alias against the grid and inject spurious oscillations
(~0.14 dex for $\sigma_{\rm UV} = 0.2$ at `dlog10m=0.1`). `CrossCorrelation`
warns at construction if the kernel spans fewer than about five grid cells.

The public API is h-free (Mpc, $M_\odot$, Mpc$^{-3}$); the calculation runs in
halomod's h-units internally and converts at the boundary.

The $\xi(r)$ spline is cached against the HOD parameter values, not the object
identity, because callers routinely mutate the model in place via
`update_parameters` inside MCMC loops.

## Example

```python
import numpy as np
from halogal import HODModel
from halogal.crosscorr import Selection, CrossCorrelation

hod = HODModel(z=8.0, eps0=0.25, Mc=10**11.98, sigma_UV=0.5, Mcut=1e10)

hosts    = Selection(-22.0, label="Euclid dropouts")
emitters = Selection(-20.0, -22.0, label="[OIII] emitters")

xc = CrossCorrelation(hod, hosts, emitters)

print(xc.number_density(hosts), xc.bias_effective(hosts))

# projected, both samples spectroscopic
wp = xc.wp(np.logspace(-2, 1, 20), pi_max=10.0)

# angular, spec-z hosts against a broad photometric sample
spec  = lambda z: np.where((z >= 6.75) & (z <= 9.05), 1.0, 0.0)
photo = lambda z: np.exp(-0.5 * ((z - 8.0) / 1.2) ** 2)
w = xc.w_theta(np.logspace(0, 2.2, 20), spec, photo, 5.0, 12.0)

# expected companions within 1 cMpc, +/- 10 cMpc along the line of sight
n_expected, n_random = xc.mean_companions(rp=1.0, pi_max=10.0)
```

## Caveats

- $R_{cs} = 0$ (independent scatter) is an assumption — see above.
- $\sigma_{\rm UV}$ is constant with halo mass in halogal. Models with
  $\sigma_{\rm UV}(M_h)$ rising toward low mass behave differently.
- $M_{\rm cut}$ directly sets $N_s^2$ and therefore the one-halo amplitude. It
  stops being a nuisance parameter for cross-correlation work and must be
  marginalised over.
- The Limber $w(\theta)$ evaluates $\xi$ at a single redshift.
- halogal's redshift-evolution defaults are fitted at $z \sim 4$–6; at $z = 8$
  several parameters (notably $M_{\rm cut}$) are extrapolations and should be
  set explicitly rather than taken from `DEFAULT_REDSHIFT_EVOLUTION`.
