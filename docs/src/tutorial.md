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
H, D = 15.0, 4.0           # borehole depth [m] and buried depth [m]

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

## 3. The high-level interface

[`ground_response`](@ref) wraps the model in an `AbstractGroundModel` object and dispatches
automatically. For a single borehole pass a `1×2` coordinate matrix:

```julia
xy = [0.0 0.0]
m = FLSModel(H, D, ks, Cs)              # fls model
g = ground_response(t, rb, xy, m)       # ≈ fls(t, rb, H, D, ks, Cs)
```

By default (`interp = true`) the single-borehole path sub-samples the model and interpolates back, so `g` is very close to (but not bit-identical to) a direct `fls` call. Pass `interp = false` for the exact evaluation at every `t`.

The advantage is that the *same call* works for a field, only the coordinate matrix changes.

## 4. A borehole field

Generate a layout with [`borefield`](@ref), then call [`ground_response`](@ref): it detects the multiple boreholes and applies spatial superposition for you.

```julia
xy = borefield(:rectangle, 3, 3, 5.0)       # 3×3 grid, 5 m spacing → 9×2 matrix
gf = ground_response(t, rb, xy, m)          # BC-II spatial superposition applied internally
```

The field response `gf` is larger than the single-borehole `g`: each borehole is warmed by its
neighbours, so the field needs a larger temperature change to exchange the same heat per metre. See
[Spatial superposition](@ref) and [Borefields](@ref) for the layouts and methods.

## 5. Choosing the boundary condition

The `bc` keyword of [`ground_response`](@ref) chooses how the load is shared among the boreholes;
each boundary condition also has its own function:

```julia
g1 = ground_response(t, rb, xy, m; bc = :I)          # = bc_I(t, rb, xy, m): same flux everywhere
g2 = ground_response(t, rb, xy, m)                   # = bc_II(t, rb, xy, m): equal mean temperature
g3 = ground_response(t, rb, xy, FLSModel(H, D, ks, Cs, 8))   # = bc_III(...): 8 segments per borehole
# g1 ≥ g2 ≥ g3 at every time
```

BC-II (the default) is the usual choice for boreholes connected in parallel. BC-III refines it by
letting the flux vary along each borehole, and needs an [`FLSModel`](@ref) with `nseg > 1`. BC-II
and BC-III are solved directly by time marching (no iteration).

!!! note "Time sampling and `interp`"
    BC-II and BC-III superpose the load history of every borehole, which requires a uniform time
    grid `t = Δt, 2Δt, …`. So that a log-spaced `t` (the usual choice) still gives a correct
    result, the `interp = true` default solves on internal uniform blocks (hours to decades) and
    interpolates back to your `t`; this also makes the cost independent of `length(t)`. Pass
    `interp = false` only when `t` is already `Δt, 2Δt, …, nΔt`. For BC-I and a single borehole,
    `interp` is a pure speed-up rather than a correctness requirement.

## 6. Adding groundwater flow

Switch to a moving model to include advection. The moving models take a Darcy velocity ``v_D`` and
the groundwater heat capacity ``C_f``. Because flow (along ``+x``) breaks radial symmetry, their
response depends on both the separation and the flow-relative angle. Through
[`ground_response`](@ref) you still just pass the coordinate array `xy`: the distances and angles are
computed from it internally by [`borefield_geometry`](@ref).

```julia
Cf, vD = 4.2e6, 1e-6
m_mils = MILSModel(rb, ks, Cs, Cf, vD)
m_mfls = MFLSModel(H, rb, D, ks, Cs, Cf, vD)

g_mils = ground_response(t, rb, xy, m_mils)
g_mfls = ground_response(t, rb, xy, m_mfls)
```

Unlike the conductive models, the moving g-functions reach a **steady state** at long times because
groundwater carries the injected heat away. See [Moving-source models](@ref).

!!! note  "Moving models need a positive velocity"
    [`mils`](@ref) and [`mfls`](@ref) divide by a Bessel function of ``v_D`` and require
    ``v_D > 0``. For near-impervious ground use a small value (e.g. `1e-12`) or use the
    corresponding non-moving model.

## 7. Working with low-level g-arrays

To build a field response from a precomputed array, for example, to reuse it across boundary
conditions, evaluate a model on the geometry matrices from [`borefield_geometry`](@ref). BC-II needs
the array on a uniform grid `t = Δt, 2Δt, …`:

```julia
tu = 2_592_000.0 .* (1:240)               # monthly steps over 20 years [s]
r = borefield_geometry(xy, rb)[1]         # nb×nb distance matrix
g3D = fls(tu, r, H, D, ks, Cs)            # nt × nb × nb g-array
g1 = bc_I(g3D)                            # field g-functions
g2 = bc_II(g3D)

# For a moving model, also pass the angle matrix:
r, θ = borefield_geometry(xy, rb)
g3D = mfls(tu, r, θ, H, rb, D, ks, Cs, Cf, vD)
g2 = bc_II(g3D)
```

## Validation scripts

The `script/` directory contains runnable scripts that exercise every model and method. Run them
from the package root:

```
julia --project=script/ -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=script/ script/script_ground_models.jl
```

| Script | What it shows |
|---|---|
| `script_ground_models.jl` | Single-borehole comparison of ILS, ICS, FLS, MILS, MFLS on a log-time axis. |
| `script_spatial_superposition.jl` | BC-I / BC-II / BC-III on a 3×3 field; BC-III convergence with segments; scaling from 1×1 to 5×5; computing time of `bc_II` / `bc_III` against the legacy solvers; `interp = true` against an hourly direct solve. |
| `script_borefield_layouts.jl` | All six layout configurations; g-functions of two fields rotated from 0 to 90° against the groundwater flow (FLS and MFLS). |
| `script_ground_response.jl` | `ground_response` dispatch across all models, with `@btime` benchmarks. |
| `script_groundwater_advection.jl` | 2-D heatmaps of FLS, MILS, MFLS at 10 yr — the asymmetric plume. |
