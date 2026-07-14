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

ks, Cs = 3.0, 2.0e6          # ground conductivity [W/m·K], heat capacity [J/m³·K]
rb     = 0.076               # borehole radius [m]

t = 3600.0 .* exp10.(range(0, log10(8760 * 25), length = 200))   # 1 h → 25 yr [s]
```

## 2. A single borehole, raw functions

The raw functions evaluate a model directly at a radius. For the finite line source you also supply
the depth ``H`` and the buried depth ``D``:

```julia
g_ils = ils(t, rb, ks, Cs)              # infinite line source
g_ics = ics(t, rb, rb, ks, Cs)          # infinite cylindrical source (rc = rb)
g_fls = fls(t, rb, 150.0, 4.0, ks, Cs)  # finite line source, H = 150 m, D = 4 m
```

At short times the ICS sits below the ILS (it accounts for the finite borehole radius); at long
times the FLS levels off to a steady state while the ILS keeps growing. See
[Line-source models](@ref) for the equations behind each.

## 3. The high-level interface

[`ground_response`](@ref) wraps the model in an `AbstractGroundModel` object and dispatches
automatically. For a single borehole pass a `1×2` coordinate matrix:

```julia
m = FLSModel(150.0, 4.0, ks, Cs)            # H = 150 m, D = 4 m
g = ground_response(t, rb, [0.0 0.0], m)    # ≈ fls(t, rb, 150.0, 4.0, ks, Cs)
```

By default (`interp = true`) the single-borehole path sub-samples the model and interpolates back,
so `g` is very close to — but not bit-identical to — a direct `fls` call. Pass `interp = false` for
the exact evaluation at every `t`.

The advantage is that the *same call* works for a field — only the coordinate matrix changes.

## 4. A borehole field

Generate a layout with [`borefield`](@ref), then call [`ground_response`](@ref): it detects the
multiple boreholes and applies spatial superposition for you.

```julia
xy = borefield(:rectangle, 3, 3, 5.0)       # 3×3 grid, 5 m spacing → 9×2 matrix
gf = ground_response(t, rb, xy, m)          # successive_flux applied internally
```

The field response `gf` is larger than the single-borehole `g`: the boreholes warm each other's
ground. See [Spatial superposition](@ref) and [Borefields](@ref) for the layouts and methods.

## 5. Choosing the superposition method

`ground_response` uses [`successive_flux`](@ref) (fast, iterative). For a reference solution, or to
validate it, call [`bloc_matrix`](@ref) (a direct linear solve). Both accept either a model or a
precomputed g-array:

```julia
g_succ = successive_flux(t, rb, xy, m)
g_bloc = bloc_matrix(t, rb, xy, m)
# g_succ ≈ g_bloc — the two methods agree to within iteration tolerance
```

!!! note "Time sampling and `interp`"
    [`successive_flux`](@ref) reproduces the load history with a convolution that assumes a
    **constant time step**. So that a log-spaced `t` (the usual choice) still gives a correct
    result, the `interp = true` default solves on an internal constant-step grid and interpolates
    back to your `t` — this also makes the cost independent of `length(t)`. Pass `interp = false`
    only when `t` is already uniformly spaced. The instantaneous methods ([`bloc_matrix`](@ref),
    [`uniform_flux`](@ref), and the single-borehole path) are spacing-independent, so for them
    `interp` is a pure speed-up rather than a correctness requirement.

## 6. Adding groundwater flow

Switch to a moving model to include advection. The moving models take a Darcy velocity ``v_D`` and
the groundwater heat capacity ``C_f``. Because flow (along ``+x``) breaks radial symmetry, their
response depends on both the separation and the flow-relative angle — but through
[`ground_response`](@ref) you still just pass the coordinate array `xy`:

```julia
Cf, vD = 4.2e6, 1e-6
m_mils = MILSModel(rb, ks, Cs, Cf, vD)
m_mfls = MFLSModel(150.0, rb, 4.0, ks, Cs, Cf, vD)

g_mils = ground_response(t, rb, xy, m_mils)
g_mfls = ground_response(t, rb, xy, m_mfls)
```

Unlike the conductive models, the moving g-functions reach a **steady state** at long times because
groundwater carries the injected heat away. See [Moving-source models](@ref).

!!! warning "Moving models need a positive velocity"
    [`mils`](@ref) and [`mfls`](@ref) divide by a Bessel function of ``v_D`` and require
    ``v_D > 0``. For near-impervious ground use a small value (e.g. `1e-12`) or use the
    corresponding non-moving model.

## 7. Working with low-level g-arrays

To build a field response from a precomputed array — for example, to reuse it across superposition
methods — evaluate a model on the geometry matrices from [`borefield_geometry`](@ref):

```julia
r3D = borefield_geometry(xy, rb)[1]         # nb×nb distance matrix
g3D = fls(t, r3D, 150.0, 4.0, ks, Cs)       # nt × nb × nb g-array
g   = successive_flux(g3D)                  # field g-function

# For a moving model, also pass the angle matrix:
r, θ = borefield_geometry(xy, rb)
g3D  = mfls(t, r, θ, 150.0, rb, 4.0, ks, Cs, Cf, vD)
g    = successive_flux(g3D)
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
| `script_spatial_superposition.jl` | Block matrix vs successive flux agreement (3×3 FLS); scaling from 1×1 to 5×5. |
| `script_borefield_layouts.jl` | A subplot figure of all six layout configurations. |
| `script_ground_response.jl` | `ground_response` dispatch across all models, with `@btime` benchmarks. |
| `script_groundwater_advection.jl` | 2-D heatmaps of FLS, MILS, MFLS at 10 yr — the asymmetric plume. |
| `script_short_term.jl` | Short-term ANN model (⚠️ not yet fully integrated). |
