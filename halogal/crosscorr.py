"""Halo-model cross-correlations between two magnitude-selected samples.

This module computes the cross-correlation between two UV-magnitude-selected
galaxy samples drawn from the same underlying HOD -- for example a sparse,
bright, spectroscopically confirmed sample and a fainter, denser sample found
in the same fields.

The implementation is direct (a "Route B" calculation): occupations come from
:class:`~halogal.model.HODModel`, while the halo mass function, halo bias,
tracer profile and halo-centre power spectrum are taken from a single halomod
``TracerHaloModel`` backbone. Using one backbone for both samples guarantees
that every ingredient shares a mass definition and an HMF -- mixing sources
here is a known way to introduce ~0.5 dex errors.

Equations
---------
For a sample :math:`X` with occupations :math:`N_c^X(M)`, :math:`N_s^X(M)`,

.. math:: \\bar n_X = \\int d\\ln M\\, \\frac{dn}{d\\ln M}
          \\left[N_c^X(M) + N_s^X(M)\\right].

The one-halo cross power spectrum is

.. math:: P^{1h}_{12}(k) = \\frac{1}{\\bar n_1 \\bar n_2}
          \\int d\\ln M\\, \\frac{dn}{d\\ln M}
          \\Big[\\langle N_c^1 N_s^2\\rangle u(k|M)
              + \\langle N_s^1 N_c^2\\rangle u(k|M)
              + \\langle N_s^1 N_s^2\\rangle u^2(k|M)\\Big],

with no central-central term: a halo hosts one central, so for disjoint
selections that central belongs to at most one sample. The pair counts allow
for correlated scatter between the two samples' occupations,

.. math:: \\langle N_c^1 N_s^2\\rangle = N_c^1 N_s^2
          + R_{cs}\\,\\sigma_c^1 \\sigma_s^2,

where :math:`\\sigma_c = \\sqrt{N_c(1-N_c)}` (Bernoulli) and
:math:`\\sigma_s = \\sqrt{N_s}` (Poisson). ``R_cs = R_sc = R_ss = 0``
(the default) means the scatter in :math:`M_{\\rm UV}` is independent between
a central and its satellites.

The two-halo term is

.. math:: P^{2h}_{12}(k) = b_1(k)\\, b_2(k)\\, P_{hh}(k), \\quad
          b_X(k) = \\frac{1}{\\bar n_X}\\int d\\ln M\\, \\frac{dn}{d\\ln M}
          b(M)\\left[N_c^X + N_s^X u(k|M)\\right].

The real-space correlation function follows by Hankel transform, and the two
projections are

.. math:: w_p(r_p) = 2\\int_0^{\\pi_{\\rm max}} d\\pi\\,
          \\xi_{12}\\!\\left(\\sqrt{r_p^2 + \\pi^2}\\right),

.. math:: w_{12}(\\theta) = \\int d\\chi\\, p_1(\\chi)\\, p_2(\\chi)
          \\int du\\, \\xi_{12}\\!\\left(\\sqrt{u^2 + \\chi^2\\theta^2}\\right).

See ``docs/crosscorrelations.md`` for the full derivation and the modelling
caveats.

Conventions that differ from halomod
------------------------------------
Two deliberate differences, both validated against halomod 2.2.3:

1. In the two-halo term the profile transform :math:`u(k|M)` multiplies only
   the **satellite** occupation. A central sits at the halo centre, so its
   transform is unity by definition. halomod instead applies ``u`` to the
   total occupation (``cross_correlations.py`` ``power_2h_cross_fnc``), which
   makes its :math:`b(k)` collapse at high :math:`k`. The two agree to
   :math:`10^{-4}` for :math:`k \\lesssim 0.3\\,h\\,{\\rm Mpc}^{-1}` and diverge
   above :math:`k \\sim 1`, where the one-halo term dominates by orders of
   magnitude, so the total power still agrees to :math:`5\\times10^{-4}`.

2. halomod's ``CrossCorrelations.corr_cross`` returns :math:`1 + \\xi`, not
   :math:`\\xi` -- it adds a spurious ``+ 1`` that its own auto-correlation
   path does not. ``CrossCorrelation.xi`` returns :math:`\\xi`. Subtract one
   from ``corr_cross`` before comparing, or large-separation agreement will
   look broken.

**Halo exclusion is not modelled.** halomod's default ``NgMatched_`` exclusion
suppresses the two-halo term below roughly the sum of two virial radii; this
module has no exclusion, so the one-to-two-halo transition is sharper than
``halogal.model.Observables`` would give for the same HOD.

Units
-----
The public API is h-free: masses in :math:`M_\\odot`, separations in Mpc,
number densities in :math:`{\\rm Mpc}^{-3}`. Internally the calculation runs
in halomod's h-units and converts at the boundary.
"""

import warnings

import numpy as np
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline

from .config import DEFAULT_COSMO_PARAMS, halomod_cosmo_kwargs

try:
    import halomod
    from halomod.tools import ExtendedSpline, hankel_transform

    HALOMOD_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised only without halomod
    HALOMOD_AVAILABLE = False


__all__ = ["Selection", "CrossCorrelation"]

# np.trapz was renamed np.trapezoid in numpy 2.0; support both.
_trapz = getattr(np, "trapezoid", None) or np.trapz


class Selection:
    """A UV absolute magnitude selection.

    Selects galaxies with ``bright < MUV < faint``. Because UV absolute
    magnitudes are negative and brighter objects are more negative, ``bright``
    is the more negative bound.

    Parameters
    ----------
    faint : float
        Faint magnitude limit (the least negative bound). Galaxies fainter
        than this are excluded.
    bright : float, optional
        Bright magnitude limit (the more negative bound). Default ``-np.inf``,
        i.e. a pure threshold sample of everything brighter than ``faint``.
    label : str, optional
        Human-readable name, used in warnings and reprs.

    Examples
    --------
    >>> hosts = Selection(faint=-22.0, label="Euclid dropouts")   # MUV < -22
    >>> emitters = Selection(faint=-20.0, bright=-22.0)           # a bin

    Notes
    -----
    For a cross-correlation the two selections should be disjoint. If they
    overlap, a galaxy appears in both samples and contributes a self-pair at
    zero separation, which is not part of the correlation function.
    """

    def __init__(self, faint, bright=-np.inf, label=None):
        if bright >= faint:
            raise ValueError(
                f"bright limit ({bright}) must be more negative than the faint "
                f"limit ({faint}); magnitudes get more negative when brighter"
            )
        self.faint = float(faint)
        self.bright = float(bright)
        self.label = label

    @property
    def is_threshold(self):
        """bool: True if this is an open-ended threshold selection."""
        return not np.isfinite(self.bright)

    def Ncen(self, hod_model, Mh):
        """Mean central occupation of this selection.

        Parameters
        ----------
        hod_model : HODModel
            Model supplying ``Ncen(Mh, MUV_thresh)``.
        Mh : ndarray
            Halo masses in M_sun.

        Returns
        -------
        ndarray
            Mean number of centrals per halo falling in the selection.
        """
        n = hod_model.Ncen(Mh, self.faint)
        if not self.is_threshold:
            n = n - hod_model.Ncen(Mh, self.bright)
        return n

    def Nsat(self, hod_model, Mh):
        """Mean satellite occupation of this selection.

        Parameters
        ----------
        hod_model : HODModel
            Model supplying ``Nsat(Mh, MUV_thresh)``.
        Mh : ndarray
            Halo masses in M_sun.

        Returns
        -------
        ndarray
            Mean number of satellites per halo falling in the selection.
        """
        n = hod_model.Nsat(Mh, self.faint)
        if not self.is_threshold:
            n = n - hod_model.Nsat(Mh, self.bright)
        return n

    def overlaps(self, other):
        """bool: whether this selection shares magnitudes with ``other``."""
        return (self.bright < other.faint) and (other.bright < self.faint)

    def __repr__(self):
        rng = (
            f"MUV < {self.faint}"
            if self.is_threshold
            else f"{self.bright} < MUV < {self.faint}"
        )
        lab = f" '{self.label}'" if self.label else ""
        return f"Selection({rng}){lab}"


class CrossCorrelation:
    """Cross-correlation of two magnitude-selected samples in one HOD.

    Parameters
    ----------
    hod_model : HODModel
        The occupation model. Its HOD parameters may be updated in place
        between calls; results are recomputed on demand. **The redshift is
        the exception**: the halomod backbone is built once at
        ``hod_model.z`` and cannot follow a later change, so mutating ``z``
        raises rather than silently pairing occupations at one redshift with
        a mass function at another. Build a new ``CrossCorrelation`` instead.
    selection_1, selection_2 : Selection
        The two samples to cross-correlate. They should be disjoint; a warning
        is issued if they overlap.
    backbone : halomod.TracerHaloModel, optional
        Existing halomod model supplying the HMF, halo bias, tracer profile
        and halo-centre power. If omitted one is built from ``hod_model.z``
        and the package default cosmology.
    R_cs, R_sc, R_ss : float, optional
        Correlation coefficients between the occupations of the two samples
        (central-satellite, satellite-central, satellite-satellite). Zero
        means independent scatter. See module docstring.
    **backbone_kwargs
        Passed to ``halomod.TracerHaloModel`` when building the backbone.
        ``dlog10m`` defaults to 0.02: halomod's own default of 0.1 undersamples
        the ``sigma_UV`` erf kernel and aliases badly for small ``sigma_UV``.

    Attributes
    ----------
    backbone : halomod.TracerHaloModel
        The shared halo-model backbone.

    Examples
    --------
    >>> from halogal import HODModel
    >>> from halogal.crosscorr import Selection, CrossCorrelation
    >>> hod = HODModel(z=8.0)
    >>> xc = CrossCorrelation(hod, Selection(-22.0), Selection(-20.0, -22.0))
    >>> rp = np.logspace(-2, 1, 20)
    >>> wp = xc.wp(rp, pi_max=10.0)
    """

    #: halomod's default mass resolution aliases the sigma_UV erf kernel.
    DEFAULT_DLOG10M = 0.02

    def __init__(
        self,
        hod_model,
        selection_1,
        selection_2,
        backbone=None,
        R_cs=0.0,
        R_sc=0.0,
        R_ss=0.0,
        **backbone_kwargs,
    ):
        if not HALOMOD_AVAILABLE:
            raise ImportError("halomod is required for CrossCorrelation")

        self.hod_model = hod_model
        self.sel1 = selection_1
        self.sel2 = selection_2
        self.R_cs = R_cs
        self.R_sc = R_sc
        self.R_ss = R_ss

        if selection_1.overlaps(selection_2):
            warnings.warn(
                f"{selection_1!r} and {selection_2!r} overlap in magnitude. "
                "Galaxies in both samples contribute self-pairs at zero "
                "separation, which the cross-correlation does not model. "
                "Use disjoint selections unless you handle this yourself.",
                stacklevel=2,
            )

        self.backbone = backbone if backbone is not None else self._build_backbone(
            **backbone_kwargs
        )

        # The backbone is fixed at this redshift; see _assert_backbone_z.
        self._backbone_z = hod_model.z

        self._h = self.backbone.cosmo.H0.value / 100.0
        self._m_h = self.backbone.m                       # M_sun/h
        self._lnm = np.log(self._m_h)
        self._Mh = self._m_h / self._h                    # M_sun
        self._dndlnm = self.backbone.dndlnm               # (h/Mpc)^3
        self._bias = self.backbone.halo_bias
        self._u = self.backbone.tracer_profile_ukm        # (k, m), normalised
        self._k = self.backbone.k                         # h/Mpc

        self._xi_cache = None
        self._xi_cache_key = None

        self._check_mass_resolution()

    # ------------------------------------------------------------------
    # construction / validation
    # ------------------------------------------------------------------
    def _build_backbone(self, **kwargs):
        """Build a TracerHaloModel sharing the package default cosmology."""
        kwargs.setdefault("dlog10m", self.DEFAULT_DLOG10M)
        kwargs.setdefault("Mmin", 5)
        kwargs.setdefault("Mmax", 16)
        kwargs.setdefault("hmf_model", "Watson")
        kwargs.setdefault("bias_model", "Tinker10")
        # Match Observables.initialize_correlation_model, whose default is
        # nonlinear rather than halomod's own "linear", so the two halogal
        # entry points describe the same universe.
        kwargs.setdefault("hc_spectrum", "nonlinear")
        # One cosmology for the whole package: this carries sigma_8 and n as
        # well as the background parameters, so the backbone cannot drift onto
        # halomod's own Planck18 defaults while Colossus uses halogal's.
        for key, value in halomod_cosmo_kwargs().items():
            kwargs.setdefault(key, value)
        return halomod.TracerHaloModel(z=self.hod_model.z, **kwargs)

    def _check_mass_resolution(self):
        """Warn if the sigma_UV erf kernel is undersampled by the mass grid.

        The central occupation is an erf of width ``sigma_UV`` in magnitude.
        Translated to halo mass via dMUV/dlogM, a small ``sigma_UV`` can be
        narrower than one mass grid cell, in which case the mass integrals
        alias against the grid and the results oscillate.
        """
        muv = self.hod_model.MUV(self._Mh)
        dlog10m = np.diff(np.log10(self._m_h)).mean()
        dmuv_dlogm = np.abs(np.gradient(muv, np.log10(self._Mh))).max()
        n_per_kernel = self.hod_model.sigma_UV / (dmuv_dlogm * dlog10m)
        if n_per_kernel < 5:
            warnings.warn(
                f"The sigma_UV={self.hod_model.sigma_UV:.2f} kernel spans only "
                f"~{n_per_kernel:.1f} mass grid cells (dlog10m={dlog10m:.3f}). "
                "Mass integrals will alias; rebuild the backbone with a "
                "smaller dlog10m.",
                stacklevel=3,
            )

    # ------------------------------------------------------------------
    # occupations and sample properties
    # ------------------------------------------------------------------
    def _assert_backbone_z(self):
        """Fail loudly if the model's redshift has drifted from the backbone.

        The mass function, halo bias and profile are snapshotted from the
        backbone at construction. If ``hod_model.z`` is changed afterwards the
        occupations move but those ingredients do not, silently combining two
        redshifts -- a ~1 dex error in the number density for a change of
        z = 8 -> 5. There is no cheap way to follow the change, so refuse.
        """
        if self.hod_model.z != self._backbone_z:
            raise ValueError(
                f"hod_model.z changed from {self._backbone_z} to "
                f"{self.hod_model.z} after this CrossCorrelation was built. "
                "The halomod backbone is fixed at the construction redshift, "
                "so results would mix occupations at one redshift with a mass "
                "function at another. Build a new CrossCorrelation instead."
            )

    def _occupations(self, selection):
        """(N_cen, N_sat) arrays on the backbone mass grid.

        Every public quantity routes through here, so the redshift guard is
        applied once, centrally.
        """
        self._assert_backbone_z()
        return (
            selection.Ncen(self.hod_model, self._Mh),
            selection.Nsat(self.hod_model, self._Mh),
        )

    def _n_bar_h(self, selection):
        """Mean number density in (h/Mpc)^3 (internal units)."""
        nc, ns = self._occupations(selection)
        return integrate.simpson(self._dndlnm * (nc + ns), x=self._lnm)

    def number_density(self, selection):
        """Mean comoving number density of a selection.

        Parameters
        ----------
        selection : Selection
            The magnitude selection.

        Returns
        -------
        float
            Number density in Mpc^-3.
        """
        return self._n_bar_h(selection) * self._h**3

    def bias_effective(self, selection):
        """Large-scale effective bias of a selection.

        Parameters
        ----------
        selection : Selection
            The magnitude selection.

        Returns
        -------
        float
            Number-weighted effective bias (the ``u -> 1`` limit of ``b_X(k)``).
        """
        nc, ns = self._occupations(selection)
        w = self._dndlnm * (nc + ns)
        return integrate.simpson(w * self._bias, x=self._lnm) / integrate.simpson(
            w, x=self._lnm
        )

    # ------------------------------------------------------------------
    # power spectra
    # ------------------------------------------------------------------
    def power_1h(self):
        """One-halo cross power spectrum.

        Returns
        -------
        k : ndarray
            Wavenumbers in h/Mpc.
        p : ndarray
            :math:`P^{1h}_{12}(k)` in (Mpc/h)^3.
        """
        nc1, ns1 = self._occupations(self.sel1)
        nc2, ns2 = self._occupations(self.sel2)

        # Bernoulli centrals, Poisson satellites.
        sc1 = np.sqrt(np.clip(nc1 * (1 - nc1), 0, None))
        sc2 = np.sqrt(np.clip(nc2 * (1 - nc2), 0, None))
        ss1, ss2 = np.sqrt(np.clip(ns1, 0, None)), np.sqrt(np.clip(ns2, 0, None))

        cs = nc1 * ns2 + self.R_cs * sc1 * ss2
        sc = ns1 * nc2 + self.R_sc * ss1 * sc2
        ss = ns1 * ns2 + self.R_ss * ss1 * ss2

        # These are pair counts and cannot be negative. For a rare selection
        # N_c << 1 while sigma_c = sqrt(N_c(1-N_c)) ~ sqrt(N_c), so the
        # correction term exceeds the mean product once |R| >~ sqrt(N_c N_s),
        # which is tiny. A modest negative R would otherwise drive the
        # one-halo term -- and hence xi and mean_companions -- negative,
        # i.e. fewer than random companions inside the host's virial radius.
        if np.any(cs < 0) or np.any(sc < 0) or np.any(ss < 0):
            warnings.warn(
                "The requested correlation coefficients drive a one-halo pair "
                "count negative over part of the mass range, which is "
                "unphysical. For a rare sample the attainable range is roughly "
                "|R| <~ sqrt(N_c N_s) << 1, far below order-unity values. "
                "Clipping the affected pair counts to zero.",
                stacklevel=2,
            )
            cs, sc, ss = (np.clip(x, 0.0, None) for x in (cs, sc, ss))

        integrand = self._dndlnm * (
            cs * self._u + sc * self._u + ss * self._u**2
        )
        p = integrate.simpson(integrand, x=self._lnm, axis=-1)
        return self._k, p / (self._n_bar_h(self.sel1) * self._n_bar_h(self.sel2))

    def _bias_k(self, selection):
        """Scale-dependent effective bias b_X(k), shape (k,)."""
        nc, ns = self._occupations(selection)
        integrand = self._dndlnm * self._bias * (nc + ns * self._u)
        return integrate.simpson(integrand, x=self._lnm, axis=-1) / self._n_bar_h(
            selection
        )

    def power_2h(self):
        """Two-halo cross power spectrum.

        Returns
        -------
        k : ndarray
            Wavenumbers in h/Mpc.
        p : ndarray
            :math:`P^{2h}_{12}(k)` in (Mpc/h)^3.
        """
        p_hh = self.backbone._power_halo_centres_fnc(self._k)
        return self._k, self._bias_k(self.sel1) * self._bias_k(self.sel2) * p_hh

    def power(self):
        """Total cross power spectrum (one-halo plus two-halo).

        Returns
        -------
        k : ndarray
            Wavenumbers in h/Mpc.
        p : ndarray
            :math:`P_{12}(k)` in (Mpc/h)^3.
        """
        k, p1 = self.power_1h()
        _, p2 = self.power_2h()
        return k, p1 + p2

    # ------------------------------------------------------------------
    # correlation function and projections
    # ------------------------------------------------------------------
    def _param_key(self):
        """Hashable snapshot of everything xi(r) depends on.

        The HOD model is mutable and callers update it in place between
        evaluations, so the cache must key on parameter values rather than
        assume immutability.
        """
        m = self.hod_model
        return (
            m.z, m.eps0, m.Mc, m.a, m.b, m.sigma_UV, m.Mcut, m.Msat, m.asat,
            m.add_dust, self.R_cs, self.R_sc, self.R_ss,
            self.sel1.faint, self.sel1.bright,
            self.sel2.faint, self.sel2.bright,
        )

    def _xi_spline_h(self, r_min=1e-3, r_max=300.0, n_r=300):
        """Spline of xi_12(r) with r in Mpc/h (internal units).

        Cached against the current parameter values, so repeated projections
        in an MCMC loop pay for one Hankel transform rather than several.
        """
        key = (self._param_key(), r_min, r_max, n_r)
        if self._xi_cache is not None and self._xi_cache_key == key:
            return self._xi_cache

        k, p = self.power()
        p_fnc = ExtendedSpline(
            k, p, lower_func="power_law", upper_func="power_law"
        )
        r_tab = np.logspace(np.log10(r_min), np.log10(r_max), n_r)
        xi = hankel_transform(p_fnc, r_tab, "r")
        spline = ExtendedSpline(
            r_tab, xi, lower_func="power_law", upper_func="power_law"
        )
        self._xi_cache, self._xi_cache_key = spline, key
        return spline

    def xi(self, r):
        """Real-space cross-correlation function.

        Parameters
        ----------
        r : array_like
            Comoving separations in Mpc.

        Returns
        -------
        ndarray
            :math:`\\xi_{12}(r)`, dimensionless.
        """
        r = np.atleast_1d(r).astype(float)
        return self._xi_spline_h()(r * self._h)

    def _xi_projected_h(self, R_h, pi_max_h, xi_spline=None, n_pi=200):
        """2 * integral of xi along the line of sight, all in Mpc/h.

        Parameters
        ----------
        R_h : ndarray
            Transverse separations in Mpc/h.
        pi_max_h : float
            Line-of-sight integration half-length in Mpc/h.
        xi_spline : callable, optional
            Precomputed xi(r) spline in Mpc/h, to avoid recomputation.
        n_pi : int, optional
            Number of line-of-sight samples.

        Returns
        -------
        ndarray
            The projected integral, in Mpc/h.
        """
        if xi_spline is None:
            xi_spline = self._xi_spline_h()
        R_h = np.atleast_1d(R_h).astype(float)

        # Per-R sinh grid. A single shared grid (e.g. sqrt-spaced) has its
        # first node at a fixed pi_1; once R falls below it the steep core of
        # xi falls between nodes and Simpson overshoots, returning a NEGATIVE
        # w_p. Scaling the grid with R resolves both the |pi| <~ R core and
        # the long tail, at the same cost.
        scale = np.clip(R_h, pi_max_h / 1e6, pi_max_h)[:, None]
        t = np.linspace(0.0, 1.0, n_pi)[None, :]
        pi = scale * np.sinh(t * np.arcsinh(pi_max_h / scale))
        rr = np.sqrt(R_h[:, None] ** 2 + pi**2)
        return 2.0 * integrate.simpson(xi_spline(rr), x=pi, axis=-1)

    def wp(self, rp, pi_max=100.0):
        """Projected cross-correlation function.

        .. math:: w_p(r_p) = 2\\int_0^{\\pi_{\\rm max}} d\\pi\\,
                  \\xi_{12}\\!\\left(\\sqrt{r_p^2 + \\pi^2}\\right)

        Parameters
        ----------
        rp : array_like
            Transverse comoving separations in Mpc.
        pi_max : float, optional
            Line-of-sight integration limit in Mpc. Default 100.

        Returns
        -------
        ndarray
            :math:`w_p(r_p)` in Mpc.
        """
        rp = np.atleast_1d(rp).astype(float)
        wp_h = self._xi_projected_h(rp * self._h, pi_max * self._h)
        return wp_h / self._h

    def w_theta(self, theta, p_z_1, p_z_2, z_min, z_max, n_z=400,
                pi_max=200.0):
        """Angular cross-correlation function via the Limber approximation.

        .. math:: w_{12}(\\theta) = \\int d\\chi\\, p_1(\\chi) p_2(\\chi)
                  \\int du\\, \\xi_{12}\\!\\left(\\sqrt{u^2+\\chi^2\\theta^2}\\right)

        The redshift distributions are normalised internally, so only their
        shapes matter. Cross-correlating a narrow spectroscopic distribution
        with a broad photometric one dilutes the signal by their overlap --
        that dilution is captured here.

        Parameters
        ----------
        theta : array_like
            Angular separations in arcseconds.
        p_z_1, p_z_2 : callable
            Redshift distributions of the two samples, ``p(z)``. Need not be
            normalised.
        z_min, z_max : float
            Redshift range over which to integrate.
        n_z : int, optional
            Number of redshift samples. Default 400. A sharp-edged (top-hat)
            p(z) is sensitive to whether grid points land on the window
            edges; at n_z=100 that shows up as a ~1.5% oscillation in
            w(theta). A smooth p(z) converges by n_z=100.
        pi_max : float, optional
            Line-of-sight integration limit in Mpc. Default 200.

        Returns
        -------
        ndarray
            :math:`w(\\theta)`, dimensionless.

        Notes
        -----
        The correlation function is evaluated at the model's own redshift
        rather than being evolved across the window, which is appropriate for
        a narrow selection but approximate for a broad one.
        """
        theta = np.atleast_1d(theta).astype(float)
        theta_rad = theta / 3600.0 * np.pi / 180.0

        z = np.linspace(z_min, z_max, n_z)
        chi = self.backbone.cosmo.comoving_distance(z).value * self._h  # Mpc/h

        # broadcast_to: a constant p(z) (e.g. `lambda z: 1.0`) is the most
        # natural way to write a flat window, but returns a scalar.
        p1 = np.broadcast_to(np.asarray(p_z_1(z), dtype=float), z.shape)
        p2 = np.broadcast_to(np.asarray(p_z_2(z), dtype=float), z.shape)
        if np.any(p1 < 0) or np.any(p2 < 0):
            raise ValueError("redshift distributions must be non-negative")

        # Normalise in chi, so that int p(chi) dchi = 1.
        p1 = p1 / _trapz(p1, chi)
        p2 = p2 / _trapz(p2, chi)

        xi_spline = self._xi_spline_h()
        w = np.zeros_like(theta_rad)
        for i, th in enumerate(theta_rad):
            proj = self._xi_projected_h(
                chi * th, pi_max * self._h, xi_spline=xi_spline
            )
            w[i] = _trapz(p1 * p2 * proj, chi)
        return w

    def mean_companions(self, rp, pi_max, selection=None):
        """Mean number of companions within a cylinder around a sample-1 object.

        This is the estimator that stacked counts around a small number of
        known hosts actually measure:

        .. math:: \\langle N \\rangle = \\bar n_2
                  \\left[\\pi r_p^2 \\cdot 2\\pi_{\\rm max}
                  + 2\\pi \\int_0^{r_p} R\\, w_p(R)\\, dR\\right]

        Parameters
        ----------
        rp : float
            Cylinder radius in comoving Mpc.
        pi_max : float
            Cylinder half-length in comoving Mpc.
        selection : Selection, optional
            The companion sample. Defaults to ``selection_2``.

        Returns
        -------
        n_expected : float
            Mean number of companions including clustering.
        n_random : float
            Mean number expected for an unclustered sample, for reference.
        """
        sel = selection if selection is not None else self.sel2
        n_bar = self.number_density(sel)

        n_random = n_bar * np.pi * rp**2 * 2.0 * pi_max
        R = np.linspace(1e-4, rp, 200)
        excess = n_bar * 2.0 * np.pi * integrate.simpson(
            R * self.wp(R, pi_max=pi_max), x=R
        )
        return n_random + excess, n_random

    def __repr__(self):
        return (
            f"CrossCorrelation(z={self.hod_model.z}, {self.sel1!r} x "
            f"{self.sel2!r})"
        )
