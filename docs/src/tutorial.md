# Tutorial

This tutorial builds a ground thermal response from a single borehole up to a borehole field, and
shows how to switch models and add groundwater flow. Every function used here is documented in the
[API reference](@ref).

## 1. Ground properties and a time vector

The models take the ground thermal conductivity ``k_s`` and the volumetric heat capacity ``C_s``;
their ratio is the diffusivity ``\alpha = k_s/C_s``. g-functions span many orders of magnitude in
time, so a **logarithmically spaced** time vector is the natural choice:

```julia
using GroundResponse

ks, Cs = 3.0, 2.0e6         # ground conductivity [W/m·K], heat capacity [J/m³·K]
rb = 0.076                  # borehole radius [m]
H, D = 150.0, 4.0          # borehole depth [m] and buried depth [m]

t = 3600.0 .* exp10.(range(0, log10(8760 * 25), length = 200))   # 1 h → 25 yr [s]
```

## 2. A single borehole, raw functions

The raw functions evaluate a model directly at a radius.

```julia
g_ils = ils(t, rb, ks, Cs)              # infinite line source
g_ics = ics(t, rb, rb, ks, Cs)          # infinite cylindrical source (rc = rb)
g_fls = fls(t, rb, H, D, ks, Cs)        # finite line source
```

At short times the ICS sits below the ILS (it accounts for the finite borehole radius); at long
times the FLS levels off to a steady state while the ILS keeps growing. See
[Line-source models](@ref) for the equations behind each.

## 3. The model types

[`gmatrix`](@ref) wraps the model in an `AbstractGroundModel` object and dispatches automatically.
For a single borehole pass a `1×2` coordinate matrix:

```julia
xy = [0.0 0.0]
m = FLSModel(H, D, ks, Cs)              # fls model
g = gmatrix(t, rb, xy, m)               # = fls(t, rb, H, D, ks, Cs)
```

The advantage is that the *same call* works for a field, only the coordinate matrix changes.

## 4. A borehole field

Generate a layout with [`borefield`](@ref), then call [`gmatrix`](@ref): it detects the multiple
boreholes and returns the pairwise response array.

```julia
xy = borefield(:rectangle, 3, 3, 5.0)       # 3×3 grid, 5 m spacing → 9×2 matrix
G = gmatrix(t, rb, xy, m)                   # nt × 9 × 9
```

`G[k, i, j]` is the response at borehole `i` to a unit flux (1 W/m) at borehole `j` at time `t[k]`.

The field response is larger than the single-borehole `g`: each borehole is warmed by its
neighbours, so the field needs a larger temperature change to exchange the same heat per metre. The
other boundary conditions (BC-II to BC-IV) share the load unequally among the boreholes and are
solved by [GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl). See
[Borefields](@ref) for the layouts.

## 5. Segments

The second method of [`gmatrix`](@ref), for an [`FLSModel`](@ref) only, splits each borehole into
`nseg` segments of equal length, numbered from the top:

```julia
Gs = gmatrix(t, rb, xy, m, 8)               # nt × 72 × 72, 8 segments per borehole
```

`Gs[k, p, q]` is the mean response over segment `p` to a unit flux on segment `q`, with the
segments ordered borehole-major (`p = (borehole − 1)⋅nseg + segment`). This is the array of the
BC-III and BC-IV superpositions of GroundHeatExchanger.jl, where the number of segments is the
`nseg` keyword of `bc_III`, `bc_IV` and `ground_response`.

## 6. Adding groundwater flow

Switch to a moving model to include advection. The moving models take a Darcy velocity ``v_D`` and
the groundwater heat capacity ``C_f``. Because flow (along ``+x``) breaks radial symmetry, their
response depends on both the separation and the flow-relative angle. Through
[`gmatrix`](@ref) you still just pass the coordinate array `xy`: the distances and angles are
computed from it internally by [`borefield_geometry`](@ref). The moving models carry their own
`rb`, which must match the `rb` passed to `gmatrix`.

```julia
Cf, vD = 4.2e6, 1e-6
m_mils = MILSModel(rb, ks, Cs, Cf, vD)
m_mfls = MFLSModel(H, rb, D, ks, Cs, Cf, vD)

G_mils = gmatrix(t, rb, xy, m_mils)
G_mfls = gmatrix(t, rb, xy, m_mfls)
```

Unlike the conductive models, the moving g-functions reach a **steady state** at long times because
groundwater carries the injected heat away. See [Moving-source models](@ref).

!!! note  "Moving models need a positive velocity"
    [`mils`](@ref) and [`mfls`](@ref) divide by a Bessel function of ``v_D`` and require
    ``v_D > 0``. For near-impervious ground use a small value (e.g. `1e-12`) or use the
    corresponding non-moving model.

## 7. Working with low-level g-arrays

[`gmatrix`](@ref) is a thin layer over the kernels: the same arrays come from evaluating a model on
the geometry matrices of [`borefield_geometry`](@ref):

```julia
r = borefield_geometry(xy, rb)[1]         # nb×nb distance matrix
G = fls(t, r, H, D, ks, Cs)               # = gmatrix(t, rb, xy, m)

# For a moving model, also pass the angle matrix:
r, θ = borefield_geometry(xy, rb)
G = mfls(t, r, θ, H, rb, D, ks, Cs, Cf, vD)   # = gmatrix(t, rb, xy, m_mfls)

# Segment-to-segment FLS: top half of a borehole on its bottom half, 5 m apart
g = fls(t, 5.0, H / 2, D, H / 2, D + H / 2, ks, Cs)
```

## Validation scripts

The `script/` directory contains runnable scripts that exercise every model. Run them from the
package root:

```
julia --project=script/ -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=script/ script/script_ground_models.jl
```

| Script | What it shows |
|---|---|
| `script_ground_models.jl` | Single-borehole comparison of ILS, ICS, FLS, MILS, MFLS on a log-time axis. |
| `script_groundwater_advection.jl` | 2-D heatmaps of FLS, MILS, MFLS at 10 yr — the asymmetric plume. |
