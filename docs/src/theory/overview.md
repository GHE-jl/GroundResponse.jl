# Overview

This page describes the role of a ground thermal response function, the conventions every model in
the package shares, and the high-level [`ground_response`](@ref) interface that dispatches over
model type and field size. The subsequent theory pages derive each model in turn.

## The ground response function

Once heat leaves the borehole wall it diffuses into the surrounding ground, whose temperature rises
slowly over months and years. The **ground thermal response function** ``g(t)`` captures this
transient: for a constant heat extraction/injection rate ``q`` [W/m] applied from ``t = 0``, the
temperature rise at a radius ``r`` from the source is

```math
\Delta T(r, t) = q\, g(r, t).
```

The single-borehole model functions return ``g`` in units of **°C·m/W**, with the
ground-conductivity normalisation (the ``1/(4\pi k_s)`` exponential-integral prefactor and its
analogues) already folded into the returned value — so no extra scaling is needed to obtain a
temperature from a load. A constant unit impulse of 1 W/m therefore yields ``\Delta T = g``.

Each model is the solution of the transient heat-conduction equation for a particular idealisation
of the heat source:

- a **line** of infinite length (ILS),
- a hollow **cylinder** of finite radius (ICS),
- a **line of finite length** that accounts for the ground surface and the borehole bottom (FLS),
- and the two **moving** variants (MILS, MFLS) that add groundwater advection.

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

The `nb×nb` pairwise-radius matrix is produced by [`borefield_radius`](@ref).

### Moving models — directional coordinates `xy`

Because groundwater flow breaks radial symmetry, [`mils`](@ref) and [`mfls`](@ref) take Cartesian
coordinates instead of a radius (flow is along the positive ``x``-axis):

| `t` | `xy` | Output | Meaning |
|---|---|---|---|
| `Real` | `AbstractVector` (`[x, y]`) | scalar | one time, one point |
| `AbstractVector` | `AbstractVector` (`[x, y]`) | `nt`-vector | time series at one point |
| `Real` | `AbstractMatrix` (`nb×2`) | `nb × nb` matrix | borefield, one time |
| `AbstractVector` | `AbstractMatrix` (`nb×2`) | `nt × nb × nb` array | full borefield g-array |

For the matrix overloads the kernel internally forms pairwise displacements; the diagonal
(self-response) uses ``[0, 0]``, which triggers the inside-borehole branch.

## The high-level interface

[`ground_response`](@ref) is the single entry point. It dispatches on the ground-model type and on
the number of boreholes:

- **Single borehole** (`size(xy, 1) == 1`): evaluates the model kernel directly at the borehole
  wall radius `rb`.
- **Multiple boreholes**: applies spatial superposition via [`successive_flux`](@ref).

```julia
m = FLSModel(150.0, 4.0, 3.0, 2.0e6)
g_single = ground_response(t, rb, [0.0 0.0], m)        # single borehole
g_field  = ground_response(t, rb, borefield(:rectangle, 3, 3, 5.0), m)  # borefield
```

This keeps user code independent of which model is used: swapping `FLSModel` for `MFLSModel`
changes only the model object, not the call.

## Extending with custom models

`AbstractGroundModel` is the extension point. Subtype it and add
[`successive_flux`](@ref) / [`bloc_matrix`](@ref) overloads to make a new model usable everywhere
`ground_response` is:

```julia
struct MyModel <: AbstractGroundModel
    ks::Float64
    Cs::Float64
end

function GroundResponse.successive_flux(t, rb, xy, m::MyModel)
    r, = borefield_radius(xy, rb)
    return successive_flux(my_gfunc(t, r, m.ks, m.Cs))
end
```

## Functions on this page

```@docs
AbstractGroundModel
ground_response
```
