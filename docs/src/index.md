# GroundResponse.jl

*Analytical ground thermal response functions (g-functions) for boreholes and borehole fields, in pure Julia.*

`GroundResponse.jl` evaluates the **ground thermal response**, or *g*-function, that links a heat
load injected into the ground to the resulting temperature rise at the borehole wall. It is the
**transient-ground layer** of a geothermal simulation stack: it provides the response of the soil
*outside* the borehole, where [`BoreholeResistance.jl`](https://github.com/GHE-jl/BoreholeResistance.jl) stops.

The package implements five analytical models, each available both as a **raw function** and as an
`AbstractGroundModel` **type** for the high-level interface:

| Model | Function | Type | Captures |
|---|---|---|---|
| Infinite line source | [`ils`](@ref) | [`ILSModel`](@ref) | radial conduction, long borehole |
| Infinite cylindrical source | [`ics`](@ref) | [`ICSModel`](@ref) | finite borehole radius |
| Finite line source | [`fls`](@ref) | [`FLSModel`](@ref) | finite depth + axial end effects |
| Moving infinite line source | [`mils`](@ref) | [`MILSModel`](@ref) | groundwater advection |
| Moving finite line source | [`mfls`](@ref) | [`MFLSModel`](@ref) | groundwater advection + finite depth |

\TODO Update this to include all the spatial superposition.
On top of the single-borehole models, two **spatial-superposition** methods —
[`successive_flux`](@ref) and [`bloc_matrix`](@ref) — assemble the response of an arbitrary
**borehole field**, and a family of [`borefield`](@ref) helpers generate common field layouts.

The package depends only on
[`SpecialFunctions.jl`](https://github.com/JuliaMath/SpecialFunctions.jl),
[`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl),
[`DSP.jl`](https://github.com/JuliaDSP/DSP.jl) and the `LinearAlgebra` standard library.

## What a g-function is

A *g*-function is the temperature response per unit heat load. The single-borehole models in this
package return ``g`` directly in units of **°C·m/W**, with the ground-conductivity normalisation
already folded in, so that for a constant load ``q`` [W/m] applied since ``t = 0`` the temperature
rise at the borehole wall is simply

```math
\Delta T_b(t) = q g(t).
```

Downstream packages convolve this response with a time-varying load to obtain the full
borehole-wall temperature history, see [Ecosystem](@ref).

## Installation

The package is not yet registered. Install it directly from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/GHE-jl/GroundResponse.jl")
```

or, in the Pkg REPL mode (press `]`):

```
pkg> add https://github.com/GHE-jl/GroundResponse.jl
```

## Quick start

```julia
using GroundResponse

ks, Cs = 3.0, 2.0e6          # ground conductivity [W/m·K], heat capacity [J/m³·K]
rb     = 0.076               # borehole radius [m]
t      = 3600.0 .* exp10.(range(0, log10(8760*25), length = 200))   # 1 h → 25 yr [s]

# Single borehole — finite line source
m = FLSModel(150.0, 4.0, ks, Cs)            # H = 150 m, D = 4 m
g = ground_response(t, rb, [0.0 0.0], m)    # g-function [°C·m/W]

# Borefield — spatial superposition applied automatically
xy = borefield(:rectangle, 3, 3, 5.0)       # 3×3 grid, 5 m spacing
gf = ground_response(t, rb, xy, m)
```

## Manual outline

- **[Tutorial](@ref)** — a worked single-borehole and borefield example, step by step.
- **Modeling theory** — the physics behind each model, with governing equations and references:
  - [Overview](@ref) — what a g-function is, the dispatch conventions and the high-level interface.
  - [Line-source models](@ref) — ILS, ICS and FLS for purely conductive ground.
  - [Moving-source models](@ref) — MILS and MFLS for ground with groundwater advection.
  - [Spatial superposition](@ref) — building a borefield response from single-borehole responses.
- **[Borefields](@ref Borefields)** — the layout generators and the pairwise-radius helper.
- **[API reference](@ref)** — the complete docstring reference for every exported symbol.
- **[References](@ref)** — the bibliography underpinning each model.

## Conventions used throughout

| Symbol | Meaning | Unit |
|---|---|---|
| ``g`` | Ground thermal response (per unit load) | °C·m/W |
| ``q`` | Heat load per unit borehole length | W/m |
| ``t`` | Time since load onset | s |
| ``r`` | Radial distance from the source | m |
| ``r_b`` | Borehole radius | m |
| ``H`` | Borehole (active) length | m |
| ``D`` | Buried depth (top of the active length) | m |
| ``k_s`` | Ground thermal conductivity | W/m·K |
| ``C_s`` | Ground volumetric heat capacity | J/m³·K |
| ``\alpha = k_s/C_s`` | Ground thermal diffusivity | m²/s |
| ``C_f`` | Groundwater volumetric heat capacity | J/m³·K |
| ``v_D`` | Darcy (groundwater) velocity | m/s |
| ``xy`` | Borehole coordinates, ``n_b \times 2`` | m |

!!! note "Output shape follows the input types"
    Every model function dispatches on the types of `t` and the spatial argument with no runtime
    branching: a scalar time and radius give a scalar, vectors give a vector or matrix, and an
    `nb×nb` radius matrix (or `nb×2` coordinate matrix) gives a full borefield g-array. See
    [Overview](@ref) for the dispatch table.

## Ecosystem

`GroundResponse.jl` is the transient-ground layer of a three-package geothermal stack:

| Package | Role |
|---|---|
| **BoreholeResistance.jl** | Water properties + borehole thermal resistances (inside the borehole). |
| **GroundResponse.jl** | Ground *g*-function models and spatial superposition (this package). |
| **GroundHeatExchanger.jl** | Simulation orchestration; depends on and re-exports both. |

Calling `using GroundHeatExchanger` re-exports every model and borefield function documented
here, and adds temporal superposition and fluid-temperature routines on top.
