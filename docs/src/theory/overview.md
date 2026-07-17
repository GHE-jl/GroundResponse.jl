# Overview

This page describes the role of a ground thermal response function, the conventions every model in
the package shares, and the high-level [`ground_response`](@ref) interface that dispatches over
model type and field size.

## The ground response function

Once heat leaves the borehole wall it diffuses into the surrounding ground, whose temperature rises
slowly over months and years. The **ground thermal response function** ``g(t)`` captures this
transient: for a constant heat extraction/injection rate ``q`` [W/m] applied from ``t = 0``, the
temperature rise at a radius ``r`` from the source is

```math
\Delta T(r, t) = q \cdot g(r, t).
```

The single-borehole model functions return ``g`` in units of **°C·m/W**, with the
ground-conductivity normalisation (the ``1/(4\pi k_s)`` exponential-integral prefactor and its
analogues) already folded into the returned value — so no extra scaling is needed to obtain a
temperature from a load. A constant unit impulse of 1 W/m therefore yields ``\Delta T = g``.

## What each model adds

| Model | Geometry | Extra physics over ILS |
|---|---|---|
| [`ils`](@ref) | infinite line | — (baseline) |
| [`ics`](@ref) | infinite cylinder, radius ``r_c`` | finite source radius (better at short times) |
| [`fls`](@ref) | finite line, depth ``H``, buried ``D`` | axial end effects / steady state at long times |
| [`mils`](@ref) | infinite line + flow | groundwater advection (asymmetric plume) |
| [`mfls`](@ref) | finite line + flow | advection **and** finite depth |

Short-time accuracy improves going ILS → ICS, long-time accuracy improves going ILS → FLS, and the
moving models break the radial symmetry of all three when groundwater flows.

## Dispatch conventions

All model functions are overloaded on the type of the time argument `t` and the spatial argument,
with **no runtime branching** — the output shape is determined entirely by the input types. This
lets the same function serve a single evaluation, a time series, a radial profile, or a full
borefield g-array.

### Conductive models — scalar radius `r`

[`ils`](@ref), [`ics`](@ref) and [`fls`](@ref) take a radial distance `r`:

| `t` | `r` | Output | Meaning |
|---|---|---|---|
| `Real` | `Real` | scalar | one time, one radius |
| `AbstractVector` | `Real` | `nt`-vector | time series at one radius |
| `Real` | `AbstractVector` | `nr`-vector | radial profile at one time |
| `AbstractVector` | `AbstractVector` | `nt × nr` matrix | time series at several radii |
| `Real` | `AbstractMatrix` (`nb×nb`) | `nb × nb` matrix | borefield, one time |
| `AbstractVector` | `AbstractMatrix` (`nb×nb`) | `nt × nb × nb` array | full borefield g-array |

The `nb×nb` pairwise-distance matrix is produced by [`borefield_geometry`](@ref).

### Moving models — distance `r` and flow-relative angle `θ`

Because groundwater flow breaks radial symmetry, [`mils`](@ref) and [`mfls`](@ref) take a second
geometry argument: the angle ``\theta`` (in **degrees**, ``\theta \in [0, 180]``) of the
source→receiver direction relative to the flow (``+x``). Downstream is ``\theta = 0`` (warmest),
upstream ``\theta = 180`` (coolest):

| `t` | `r`, `θ` | Output | Meaning |
|---|---|---|---|
| `Real` | `Real`, `Real` | scalar | one time, one point |
| `AbstractVector` | `Real`, `Real` | `nt`-vector | time series at one point |
| `Real` | `AbstractMatrix`, `AbstractMatrix` (`nb×nb`) | `nb × nb` matrix | borefield, one time |
| `AbstractVector` | `AbstractMatrix`, `AbstractMatrix` (`nb×nb`) | `nt × nb × nb` array | full borefield g-array |

The `r` and `θ` matrices come from [`borefield_geometry`](@ref). Its diagonal self entry
``(r_b, 0)`` flows through the kernel like any other pair and yields the borehole-wall self response,
so the diagonal is not special-cased.

## The high-level interface

[`ground_response`](@ref) is the single entry point. It dispatches on the ground-model type and on
the number of boreholes:

- **Single borehole** (`size(xy, 1) == 1`): evaluates the model kernel at the borehole wall radius
  `rb`.
- **Multiple boreholes**: applies spatial superposition via [`successive_flux`](@ref) (BC-II). \TODO or any other spatial superposition. Modify this to include the variety of methods. In fact, this section is unclear. The `bc` are defined after. Reorganize this so that it is more linear.
- **`FLSModel` with `nseg > 1`** (a field): applies the BC-III segment superposition, defaulting to
  the time-marching solver [`segment_response_marching`](@ref) (the block [`segment_response`](@ref)
  is reachable with `solver = :block`). A single borehole always evaluates the whole-borehole FLS
  kernel directly, regardless of `nseg`. See [Spatial superposition](@ref) for the boundary-condition
  hierarchy (BC-I → BC-II → BC-III).

The `bc` and `solver` keywords override these defaults (`bc = :I | :II | :III`, `solver` picking the
backend), and `interp` (default `true`) controls the internal constant-step sub-sampling, a correctness requirement for the temporal solvers on non-uniform `t`, a performance option elsewhere.
See [Spatial superposition](@ref).

```julia
m = FLSModel(150.0, 4.0, 3.0, 2.0e6)
g_single = ground_response(t, rb, [0.0 0.0], m)        # single borehole
g_field = ground_response(t, rb, borefield(:rectangle, 3, 3, 5.0), m)  # borefield (BC-II)
g_bc3 = ground_response(t, rb, borefield(:rectangle, 3, 3, 5.0), FLSModel(150.0, 4.0, 3.0, 2.0e6, 8))  # BC-III
```

This keeps user code independent of which model is used: swapping `FLSModel` for `MFLSModel`
changes only the model object, not the call.

## Extending with custom models

`AbstractGroundModel` is the extension point (in model.response.jl \TODO add the link). Subtype it and add a single `_borehole_response` method
returning the pairwise `nt × nb × nb` response array, every backend (`uniform_flux`,
`successive_flux`, `bloc_matrix`, `segment_response`) and the `ground_response` interface (including
`interp` sub-sampling) then work with no further overloads:

```julia
struct MyModel <: AbstractGroundModel
    ks::Float64
    Cs::Float64
end

function GroundResponse._borehole_response(t, rb, xy, m::MyModel)
    r = borefield_geometry(xy, rb)[1]          # nb×nb distance matrix
    return my_gfunc(t, r, m.ks, m.Cs)          # nt × nb × nb
end
```

## Functions on this page

```@docs
AbstractGroundModel
ground_response
```
