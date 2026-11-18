# Overview

This page describes the role of a ground thermal response function, the conventions every model in
the package shares, and the [`gmatrix`](@ref) response arrays of a borefield.

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

## Components of each model

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
with **no runtime branching** (the output shape is determined entirely by the input types). This
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

## Borefield response arrays

[`gmatrix`](@ref) turns a ground model into the response array of a borefield, the ``G`` of the
block-matrix formulation of Dusseault et al. (2018). It takes the time vector, the borehole radius,
the borehole coordinates and a ground model:

```julia
gmatrix(t, rb, xy, m)
gmatrix(t, rb, xy, m::FLSModel, nseg)
```

It dispatches on the number of boreholes:

- **Single borehole** (`size(xy, 1) == 1`): the borehole self response at the wall radius `rb`,
  a vector over `t`.
- **Multiple boreholes**: the pairwise response array, sized ``n_t \times n_b \times n_b``, whose
  entry `[k, i, j]` is the response at borehole ``i`` to a unit flux (1 W/m) at borehole ``j``:

```math
g_{ij}(t) = g\bigl(r_{ij}, \theta_{ij}, t\bigr), \qquad
r_{ij} = \begin{cases} \lVert \mathbf{x}_i - \mathbf{x}_j \rVert & i \ne j \\ r_b & i = j, \end{cases}
```

where ``\theta_{ij}`` is the angle of the source→receiver direction relative to the flow (``+x``);
for the isotropic models the ``\theta`` dependence drops out. Both matrices come from
[`borefield_geometry`](@ref).

The second method, for an [`FLSModel`](@ref) only, splits each borehole into `nseg` segments of
equal length and returns the ``n_t \times n_S \times n_S`` array between the ``n_S = n_b \cdot
n_{seg}`` segments, built from the segment-to-segment [`fls`](@ref) kernel. The arrays are evaluated
directly at the requested `t`.

```julia
m = FLSModel(150.0, 4.0, 3.0, 2.0e6)
xy = borefield(:rectangle, 3, 3, 5.0)
g_single = gmatrix(t, rb, [0.0 0.0], m)      # single borehole, nt-vector
G = gmatrix(t, rb, xy, m)                    # nt × 9 × 9
Gs = gmatrix(t, rb, xy, m, 8)                # nt × 72 × 72, 8 segments per borehole
```

This keeps user code independent of which model is used: swapping `FLSModel` for `MFLSModel`
changes only the model object, not the call.

The spatial superposition of these arrays into the borefield g-function under the boundary
conditions BC-I to BC-IV, with the time sub-sampling it needs, is done by
[GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl) (`ground_response`,
`bc_I` to `bc_IV`). This package only provide the calculation of analytical models, not the full borefield g-function.

## Extending with custom models

`AbstractGroundModel` is the extension point, defined in
[`src/model_response.jl`](https://github.com/GHE-jl/GroundResponse.jl/blob/main/src/model_response.jl).
Subtype it and add a single [`gmatrix`](@ref) method returning the pairwise `nt × nb × nb`
response array: the spatial superposition of GroundHeatExchanger.jl then works with no further
overloads:

```julia
struct MyModel <: AbstractGroundModel
    ks::Float64
    Cs::Float64
end

function GroundResponse.gmatrix(t, rb, xy, m::MyModel)
    r = borefield_geometry(xy, rb)[1]          # nb×nb distance matrix
    return my_gfunc(t, r, m.ks, m.Cs)          # nt × nb × nb
end
```

## Functions on this page

```@docs
AbstractGroundModel
gmatrix
```
