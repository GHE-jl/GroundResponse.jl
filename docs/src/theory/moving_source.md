# Moving-source models

When groundwater flows through the ground, heat is carried away by **advection** as well as
conduction. The plume around the borehole becomes asymmetric (warmer downstream, cooler upstream)
and, crucially, the long-time response no longer grows without bound: moving water continuously
removes the injected heat, so the g-function reaches a steady state. The two moving models add this
physics to the infinite (MILS) and finite line sources (MFLS).

Both models assume a **uniform Darcy velocity** ``v_D`` directed along the positive ``x``-axis. The
relevant transport speed is the *thermal* front velocity

```math
v_T = \frac{v_D\, C_f}{C_s},
```

which rescales the groundwater Darcy velocity by the ratio of the groundwater to the bulk-ground
volumetric heat capacity. Because both kernels divide by a Bessel function of ``v_T``, they
require ``v_D > 0``; for near-impervious ground use a small value (e.g. ``10^{-12}``) or fall back
to the corresponding non-moving model (infinite and finite line source).

## Direction dependence

The response of a moving source depends on the direction to the flow as well as on the distance.
Besides the separation ``r``, the moving models therefore take the **flow-relative angle**
``\theta`` (in degrees, ``0 \le \theta \le 180``) of the source→receiver direction, measured from the
flow (``+x``). The receiver lies at

```math
x = r\cos\theta, \qquad y = r\sin\theta,
```

and only ``r`` and the downstream coordinate ``x`` enter the kernels: ``\theta = 0`` is straight
downstream (warmest), ``\theta = 180`` straight upstream (coolest), and ``\pm\theta`` give the same
response. Both arguments are scalars for a single point, or matching ``n_b \times n_b`` matrices for a
borefield, built from the coordinates by [`borefield_geometry`](@ref) (see [Overview](@ref)).

Each kernel switches between two forms depending on whether the receiver lies inside or outside the
borehole cylinder:

- **``r \le r_b`` (at/inside the wall)** — a circumferential-average form, so the self-response of a
  borehole does not depend on an arbitrary azimuth (``\theta`` is ignored);
- **``r > r_b`` (outside)** — the full direction-dependent form, used for borehole-to-borehole
  interactions in a field.

This is what makes the matrix overloads behave correctly: the diagonal of
[`borefield_geometry`](@ref) (``r = r_b``, ``\theta = 0``) evaluates the self-response with the
inside branch, while the off-diagonal entries use the directional branch with the true separation
and angle between boreholes.

## Moving infinite line source (MILS)

The MILS extends the ILS with advection (Pasquier & Lamarche, 2022). For a line source of unit
strength in a uniform flow along ``+x``, the response at ``(x, y)`` is

```math
g_{\text{MILS}}(x, y, t) = \frac{1}{4\pi k_s}\,
\exp\!\left(\frac{x\,v_T}{2\alpha}\right)
\int_{r^2/4\alpha t}^{\infty} \exp\!\left(-\psi - \frac{r^2 v_T^2}{16\alpha^2\psi}\right)
\frac{\mathrm{d}\psi}{\psi}.
```

Only the exponential prefactor depends on the direction. Its average around the circle of radius
``r`` is ``I_0(r\,v_T/2\alpha)``, with ``I_0`` the modified Bessel function of the first kind, so the
response factors into the **circumferential-average** response ``\bar g(r, t)`` and a direction
weight:

```math
g_{\text{MILS}}(x, y, t) = \bar g(r, t)\,
\frac{\exp\!\bigl(x\,v_T / 2\alpha\bigr)}{I_0\!\bigl(r\,v_T / 2\alpha\bigr)}.
```

The exponential tilts the plume downstream; dividing by ``I_0`` keeps the azimuthal average equal to
``\bar g``. With the dimensionless groups

```math
b = \left(\frac{r\, v_T}{4 \alpha}\right)^2 = \left(\frac{r\, v_D\, C_f}{4 k_s}\right)^2, \qquad
\tau = \frac{4\alpha t}{r^2} = 4\,\mathrm{Fo},
```

the average response is a leaky-aquifer (Hantush) well function scaled by ``I_0(2\sqrt b)``:

```math
\bar g(r, t) = \frac{I_0\bigl(2\sqrt{b}\bigr)}{4\pi k_s}\, W\!\left(\frac{1}{\tau}, b\right),
\qquad
W\!\left(\frac{1}{\tau}, b\right) = \int_{1/\tau}^{\infty} \exp\!\left(-\psi - \frac{b}{\psi}\right)
\frac{\mathrm{d}\psi}{\psi}.
```

[`mils`](@ref) resolves ``\bar g`` in one of three ways:

- **``b > 1``** (fast flow or a distant receiver): the integral ``W`` is evaluated directly by
  Gauss–Kronrod quadrature (relative tolerance ``10^{-8}``), since the series below lose accuracy
  for large ``b``;
- **``b \le 1`` and ``\tau \le 1/b``** (early time): the convergent series of Eq. 20 of Pasquier &
  Lamarche (2022), built on ``E_1(1/\tau)\,I_0(2\sqrt b)``;
- **``b \le 1`` and ``\tau > 1/b``** (late time): the series of Eq. 25, built on the steady-state
  term ``2K_0(2\sqrt b)`` minus ``E_1(b\tau)\,I_0(2\sqrt b)``, with ``K_0`` the modified Bessel
  function of the second kind and ``E_1`` the exponential integral.

The series are truncated at five terms. Inside the borehole (``r \le r_b``), ``\bar g`` is evaluated
at the wall ``r_b`` and scaled by ``I_0(r\,v_T/2\alpha) / I_0(r_b\,v_T/2\alpha)``, which equals 1 at
the wall itself.

## Moving finite line source (MFLS)

The MFLS combines the finite-depth treatment of the FLS with the advection of the MILS (Guo et al.,
2020). Its integrand is the FLS depth kernel ``\Phi(s)`` (see [Line-source models](@ref)) multiplied
by an advection factor, with the thermal front velocity ``v_T`` defined above:

```math
g_{\text{MFLS}}(r, \theta, t) = \frac{\mathcal I}{4\pi k_s}\int_{1/\sqrt{4\alpha t}}^{\infty}
\exp\!\left(-\frac{v_T^2}{16\alpha^2 s^2} - r^2 s^2\right)\frac{\Phi(s)}{H s^2}\, \mathrm{d}s,
```

where the direction factor is

```math
\mathcal I =
\begin{cases}
I_0\!\bigl(r\,v_T / 2\alpha\bigr) & r \le r_b \quad\text{(circumferential average, Eq. 13)} \\[4pt]
\exp\!\bigl(r\,v_T\cos\theta / 2\alpha\bigr) & r > r_b \quad\text{(directional, Eq. 10)}
\end{cases}
```

and ``\theta`` is the flow-relative angle (see [Direction dependence](@ref)). Inside the borehole
the integral is evaluated at the wall ``r_b``. The depth kernel
``\Phi(s)`` reuses the same integrated error function ``\mathrm{ierf}`` as the FLS, so the MFLS
reduces to the FLS as ``v_D \to 0`` and to the MILS as ``H \to \infty``. [`mfls`](@ref) evaluates
the integral with Gauss–Kronrod quadrature.

## Models on this page

```@docs
MILSModel
mils
MFLSModel
mfls
```
