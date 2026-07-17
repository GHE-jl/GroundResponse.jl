# Line-source models

The three conductive models (infinite line source (ILS), infinite cylindrical source (ICS) and
finite line source (FLS)) solve the transient heat-conduction equation in a homogeneous ground
with **no groundwater flow**. They differ only in how faithfully they represent the geometry of the
borehole: a zero-radius infinite line, a finite-radius infinite cylinder, or a finite-length line.
All return the borehole-wall response normalised to a 1 W/m impulse, in °C·m/W.

Throughout, ``\alpha = k_s / C_s`` is the ground thermal diffusivity and ``r`` is the radial
distance from the source (taken as the borehole radius ``r_b`` for the self-response).

## Infinite line source (ILS)

The simplest model treats the borehole as an infinitely long line of zero radius. The temperature
response is the classical exponential-integral solution (Ingersol, 1948):

```math
g_{\text{ILS}}(r, t) = \frac{1}{4\pi k_s} E_1\!\left(\frac{r^2}{4\alpha t}\right),
```

where ``E_1`` is the exponential integral (`expinti` from `SpecialFunctions.jl` supplies
``Ei``). The single dimensionless group is the **Fourier number** ``Fo =
\alpha t / r^2``: the response depends on ``r`` and ``t`` only through this ratio.

The ILS is accurate at intermediate but overestimates the response at **short times**
(it ignores the finite borehole radius) and at **long times** for real boreholes (it ignores the
finite depth). It is implemented in [`ils`](@ref) and is also the large-radius limit of the ICS.

## Infinite cylindrical source (ICS)

The ICS replaces the line with a hollow cylinder of radius ``r_c``, removing the unphysical
short-time singularity of the ILS. The Carslaw & Jaeger (1959) solution is an integral over Bessel
functions:

```math
g_{\text{ICS}}(r, t) = \frac{1}{\pi^2 k_s}\int_0^\infty
\frac{\bigl(e^{-s^2\tilde t} - 1\bigr)\bigl(J_0(\tilde r s)\,Y_1(s) - Y_0(\tilde r s)\,J_1(s)\bigr)}
     {s^2\bigl(J_1^2(s) + Y_1^2(s)\bigr)}\, \mathrm{d}s,
```

with the dimensionless radius ``\tilde r = r / r_c`` and dimensionless time
``\tilde t = \alpha t / r_c^2``. ``J_n`` and ``Y_n`` are the Bessel functions of the first and
second kind.

The implementation in [`ics`](@ref) integrates this kernel with adaptive Gauss–Kronrod quadrature
([`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl)) and has a numerical safeguard:

- for ``\tilde r > 20`` the cylinder is indistinguishable from a line, so the function returns the
  ILS result directly;

## Finite line source (FLS)

Real boreholes have a finite active length ``H`` buried a distance ``D`` below the surface. The FLS
accounts for the ground surface (a constant-temperature boundary, modelled by an image source) and
for axial heat spreading from the borehole ends, so that the response reaches a **steady state** at
long times instead of growing without bound. The Claesson & Javed (2011) form integrates over the
borehole depth:

```math
g_{\text{FLS}}(r, t) = \frac{1}{4\pi k_s}\int_{1/\sqrt{4\alpha t}}^{\infty}
\frac{e^{-r^2 s^2}}{H s^2} \Phi(H,D,s) \mathrm{d}s,
```

with the depth integral collapsed into

```math
\Phi(H,D,s) = 2\mathrm{ierf}(Hs) + 2\mathrm{ierf}(Hs + 2Ds)
        - \mathrm{ierf}(2Hs + 2Ds) - \mathrm{ierf}(2Ds),
```

and the integrated error function

```math
\mathrm{ierf}(x) = x\mathrm{erf}(x) - \frac{1}{\sqrt{\pi}}\bigl(1 - e^{-x^2}\bigr).
```

The four ``\mathrm{ierf}`` terms encode the real source over ``[D, D+H]`` and its mirror image
across the surface. The lower integration limit ``1/\sqrt{4\alpha t}`` carries the time dependence,
so larger ``t`` integrates closer to the singularity at ``s = 0`` and yields a larger response.
[`fls`](@ref) evaluates this with Gauss–Kronrod quadrature.

The FLS is the recommended default for most simulations: it captures both the short-time behaviour
reasonably and the long-time steady state that the ILS misses.

## Models on this page

```@docs
ILSModel
ils
ICSModel
ics
FLSModel
fls
```
