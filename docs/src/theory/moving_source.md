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

\TODO: update this since models have now radius and angle as inputs.
Unlike the conductive models, the moving models take **Cartesian coordinates** ``[x, y]`` rather
than a radius, because the response depends on the angle to the flow direction. Each kernel switches
between two forms depending on whether the evaluation point lies inside or outside the borehole
cylinder ``r = \sqrt{x^2 + y^2}``:

- **``r \le r_b`` (at/inside the wall)** — a circumferential-average form, so the self-response of a
  borehole does not depend on an arbitrary azimuth;
- **``r > r_b`` (outside)** — the full direction-dependent form, used for borehole-to-borehole
  interactions in a field.

This is what makes the matrix overloads behave correctly: the diagonal (self-response) evaluates the
inside branch, while off-diagonal entries use the directional branch with the true displacement
between boreholes.

## Moving infinite line source (MILS)

The MILS extends the ILS with advection (Pasquier & Lamarche, 2022). The directional response
outside the borehole is

```math
g_{\text{MILS}}(x, y, t) = \bar g(r, t)\,
\frac{\exp\!\bigl(x\,v_T / 2\alpha\bigr)}{I_0\!\bigl(r\,v_T / 2\alpha\bigr)},
```

where ``\bar g(r, t)`` is the **circumferential-average** response and ``I_0`` is the modified
Bessel function of the first kind. The exponential factor tilts the plume downstream; dividing by
``I_0`` renormalises it so that the azimuthal average is recovered.

The average kernel ``\bar g`` is evaluated from convergent series in the dimensionless groups

```math
b = \left(\frac{r\, v_D\, C_f}{4 k_s}\right)^2, \qquad
\tau = \frac{4\alpha t}{r^2} = 4\,\mathrm{Fo},
```

using the early-time expansion (Eq. 20 of Pasquier & Lamarche) when ``\tau \le 1/b`` and the
late-time expansion (Eq. 25) otherwise, with ``I_0``, ``K_0`` and the exponential integral
supplying the closed-form terms. Inside the borehole the wall response is scaled by the ratio
``I_0(r\,v_T/2\alpha) / I_0(r_b\,v_T/2\alpha)``. See [`mils`](@ref).

\TODO add in this section the resoution of \bar{g} from Pasquier et Marcotte (the integral in the equations 1 and 3 of the paper).

## Moving finite line source (MFLS)

The MFLS combines the finite-depth treatment of the FLS with the advection of the MILS (Guo et al.,
2020). Its integrand is the FLS depth kernel ``\Phi(s)`` (see [Line-source models](@ref)) multiplied
by an advection factor:
\TODO change the U with the v_T defined at the top of the file, since it is repeated instead.
```math
g_{\text{MFLS}}(x, y, t) = \frac{\mathcal I}{4\pi k_s}\int_{1/\sqrt{4\alpha t}}^{\infty}
\exp\!\left(-\frac{U^2}{16\alpha^2 s^2} - r^2 s^2\right)\frac{\Phi(s)}{H s^2}\, \mathrm{d}s,
\qquad U = \frac{v_D\, C_f}{C_s},
```

where the direction factor is

```math
\mathcal I =
\begin{cases}
I_0\!\bigl(r\,U / 2\alpha\bigr) & r \le r_b \quad\text{(circumferential average, Eq. 13)} \\[4pt]
\exp\!\bigl(r\,U\cos\theta / 2\alpha\bigr) & r > r_b \quad\text{(directional, Eq. 10)}
\end{cases}
```

and ``\theta = \mathrm{atan}(y, x)`` is the angle to the flow direction. The depth kernel
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
