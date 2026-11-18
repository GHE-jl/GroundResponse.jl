# GroundResponse.jl

[![CI](https://github.com/GHE-jl/GroundResponse.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/GHE-jl/GroundResponse.jl/actions/workflows/CI.yml)
[![Docs: dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GHE-jl.github.io/GroundResponse.jl/dev)
[![Docs: stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://GHE-jl.github.io/GroundResponse.jl/stable)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPLv3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

A Julia package providing analytical ground thermal response functions (g-functions) for single
boreholes and the pairwise response arrays of borehole fields. GroundResponse.jl is designed as a
computational backbone: the raw model kernels (`ils`, `fls`, …) evaluate uninterpolated responses
at arbitrary time steps and radii, leaving the spatial superposition of a borefield and
system-level modelling to downstream packages such as
[GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl).

## Models

Five analytical ground response models are implemented, each as a raw function and as an
`AbstractGroundModel` subtype for `gmatrix`.

| Model | Function | Struct | Parameters |
|-------|----------|--------|------------|
| Infinite line source (Ingersol, 1948) | `ils` | `ILSModel` | `ks`, `Cs` |
| Infinite cylindrical source (Carslaw & Jaeger, 1959) | `ics` | `ICSModel` | `ks`, `Cs` |
| Finite line source (Claesson & Javed, 2011) | `fls` | `FLSModel` | `H`, `D`, `ks`, `Cs` |
| Moving infinite line source (Pasquier & Lamarche, 2022) | `mils` | `MILSModel` | `rb`, `ks`, `Cs`, `Cf`, `vD` |
| Moving finite line source (Guo et al., 2020) | `mfls` | `MFLSModel` | `H`, `rb`, `D`, `ks`, `Cs`, `Cf`, `vD` |

**Parameter legend** — `H`: borehole depth [m], `D`: buried depth [m], `rb`: borehole radius [m],
`ks`: ground thermal conductivity [W/mK], `Cs`: ground volumetric heat capacity [J/m³K], `Cf`:
groundwater volumetric heat capacity [J/m³K], `vD`: Darcy velocity [m/s].

The **MILS** and **MFLS** models are direction-dependent. Each implements a single kernel that
switches between the circumferential-average form (inside/at the borehole wall, `r ≤ rb`) and the
directional form (outside, `r > rb`). Groundwater flow is assumed along the positive x-axis.
These models require `vD > 0`; use a small value (e.g. `1e-12`) for near-impervious conditions,
or use the corresponding non-moving model instead.

## Multiple Dispatch

All functions are overloaded on the type of `t` and the spatial argument, with no runtime
branching. The output shape is determined entirely by the input types.

### ILS, ICS, FLS — distance `r`

| `t` | `r` | Output | Notes |
|-----|-----|--------|-------|
| `Real` | `Real` | scalar | Single time step, single radius |
| `AbstractVector` | `Real` | `nt`-vector | Time series at one radius |
| `Real` | `AbstractVector` | `nr`-vector | Radial profile at one time step |
| `AbstractVector` | `AbstractVector` | `nt × nr` matrix | Time series at multiple radii |
| `Real` | `AbstractMatrix` (nb×nb) | `nb × nb` matrix | Single time step, borefield radius matrix |
| `AbstractVector` | `AbstractMatrix` (nb×nb) | `nt × nb × nb` 3D array | Full borefield g-matrix |

`fls` also has a segment-to-segment method `fls(t, r, H1, D1, H2, D2, ks, Cs)` (Cimmino & Bernier,
2014) for a scalar `r` and a scalar or vector `t`: the mean response over a receiving segment
`(H2, D2)` to a unit flux on an emitting segment `(H1, D1)`.

### MILS, MFLS — distance `r` **and** flow-relative angle `θ`

The moving models are direction-dependent, so they dispatch on geometry exactly like the isotropic
models but take a second geometry argument: the angle `θ` (in **degrees**, `θ ∈ [0, 180]`) of the
source→receiver direction relative to the flow (+x). Downstream is `θ = 0` (warmest), upstream
`θ = 180` (coolest). Signatures are `mils(t, r, θ, rb, ks, Cs, Cf, vD)` and
`mfls(t, r, θ, H, rb, D, ks, Cs, Cf, vD)`.

| `t` | `r`, `θ` | Output | Notes |
|-----|----------|--------|-------|
| `Real` | `Real`, `Real` | scalar | Single time step, single point |
| `AbstractVector` | `Real`, `Real` | `nt`-vector | Time series at one point |
| `Real` | `AbstractMatrix`, `AbstractMatrix` (nb×nb) | `nb × nb` matrix | Single time step, full borefield |
| `AbstractVector` | `AbstractMatrix`, `AbstractMatrix` (nb×nb) | `nt × nb × nb` 3D array | Full borefield g-matrix |

The `nb×nb` distance matrix `r` (diagonal = `rb`) and angle matrix `θ` (diagonal = `0`) are both
produced by `borefield_geometry(xy, rb)`. The diagonal self entry `(rb, 0)` flows through the
kernel like any other pair and yields the borehole-wall self response, so no case is special-cased.

## Response Arrays

`gmatrix` turns a ground model into the response array of a borefield, the `G` of the block-matrix
formulation (Dusseault et al., 2018):

```julia
g = gmatrix(t, rb, xy, m::AbstractGroundModel)
g = gmatrix(t, rb, xy, m::FLSModel, nseg)
```

- **Single borehole** (`size(xy, 1) == 1`): the borehole self response at the wall (`r = rb`),
  returned as a vector over `t` (scalar for a scalar `t`).
- **Multiple boreholes**: `g[k, i, j]` is the response at borehole `i` to a unit flux (1 W/m) at
  borehole `j` at time `t[k]`, sized `nt × nb × nb` (`nb × nb` for a scalar `t`).
- **Segments** (second method, `FLSModel` only): each borehole is split into `nseg` segments of
  equal length, numbered from the top. `g[k, p, q]` is the mean response over segment `p` to a unit
  flux on segment `q`, sized `nt × nS × nS` with `nS = nb⋅nseg` and the segments ordered
  borehole-major (`p = (borehole − 1)⋅nseg + segment`). This is the array of the BC-III and BC-IV
  superpositions (Cimmino & Bernier, 2014; Cimmino, 2015).

The arrays are evaluated directly at the requested `t`. The moving models carry their own `rb`,
which must match the `rb` passed to `gmatrix`.

```julia
using GroundResponse

ks, Cs = 3.0, 2.0e6          # Ground thermal conductivity [W/mK], heat capacity [J/m³K]
rb     = 0.076               # Borehole radius [m]
t      = 3600.0 .* exp10.(range(0, log10(8760*25), length=200))  # 1 h → 25 yr

# Single borehole — FLS model
m = FLSModel(150.0, 4.0, ks, Cs)           # H=150 m, D=4 m
g = gmatrix(t, rb, [0.0 0.0], m)           # g-function [°Cm/W], = fls(t, rb, 150.0, 4.0, ks, Cs)

# Borefield — same model
xy = borefield(:rectangle, 3, 3, 5.0)      # 3×3 grid, 5 m spacing [m]
G  = gmatrix(t, rb, xy, m)                 # nt × 9 × 9
Gs = gmatrix(t, rb, xy, m, 8)              # nt × 72 × 72, 8 segments per borehole

# Moving model (with groundwater advection)
Cf, vD = 4.2e6, 1e-6
m_mfls = MFLSModel(150.0, rb, 4.0, ks, Cs, Cf, vD)
G = gmatrix(t, rb, xy, m_mfls)
```

## Spatial Superposition

The borefield g-function under the boundary conditions BC-I to BC-IV (`bc_I`, `bc_II`, `bc_III`,
`bc_IV`), the high-level `ground_response` and the time sub-sampling they need are in
[GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl) since GroundResponse 3.0.
They consume the arrays of `gmatrix`.

## Borefield Layout Functions

`borefield(shape, args...)` is the unified entry point that forwards arguments to the
corresponding layout function. All layout functions are also exported individually.

| Call | Layout | Boreholes |
|:-----|:-------|:----------|
| `borefield(:rectangle, nx, ny, B)` | rectangular `nx×ny` grid, uniform spacing `B` [m] | `nx × ny` |
| `borefield(:rectangle, nx, ny, Bx, By)` | rectangular grid with independent x/y spacing | `nx × ny` |
| `borefield(:line, n, B)` | single row along the x-axis, spacing `B` [m] | `n` |
| `borefield(:circle, nb, R)` | boreholes evenly distributed on a circle of radius `R` [m] | `nb` |
| `borefield(:L, n1, n2, B)` | L-shape: `n1` along x + `n2` along y, shared corner at origin | `n1+n2-1` |
| `borefield(:L, n1, n2, B1, B2)` | L-shape with independent arm spacings | `n1+n2-1` |
| `borefield(:U, nx, ny, B)` | U-shape: `nx` across the base + `ny` up each side | `nx+2(ny-1)` |
| `borefield(:U, nx, ny, Bx, By)` | U-shape with independent spacings | `nx+2(ny-1)` |
| `borefield(:open_rectangle, nx, ny, B)` | hollow rectangle, boreholes on perimeter only | `2(nx+ny-2)` |
| `borefield(:open_rectangle, nx, ny, Bx, By)` | hollow rectangle with independent spacings | `2(nx+ny-2)` |

All functions return an `nb×2` matrix of borehole coordinates `[x y]`. The `angle` keyword [°]
rotates the field counter-clockwise about the origin, e.g. `borefield(:rectangle, 3, 3, 5.0;
angle = 30)`, to orient it relative to the groundwater flow (+x) of the moving models or to the
site for a plot.

`borefield_geometry(xy, rb)` returns `(r, θ)`: the pairwise distance matrix (diagonal `rb`) and the
flow-relative angle matrix in degrees (diagonal `0`). These two matrices are what the models
consume. If you want to inspect a layout's geometric redundancy — the distinct `(r, θ)` combinations
and how many borehole pairs share each — you can derive it directly from `r` and `θ` (see the
`borefield_geometry` docstring for a copy-paste snippet).

## Extending with Custom Models

Subtype `AbstractGroundModel` and add a single `gmatrix` method returning the pairwise
`nt × nb × nb` response array. The spatial superposition of GroundHeatExchanger.jl (`bc_I`,
`bc_II` and `ground_response`) then works automatically:

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

No changes to the core package are needed. (For a direction-dependent model, key the geometry on
`(r, θ)` from `borefield_geometry` as `MILSModel`/`MFLSModel` do.)

## Scripts

Run from the package root with `julia --project=script/ script/<name>.jl`.
First-time setup:
```
julia --project=script/ -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```

| Script | Description |
|--------|-------------|
| `script_ground_models.jl` | Single-borehole g-function comparison: ILS, ICS, FLS, MILS, MFLS on a log-spaced time axis |
| `script_groundwater_advection.jl` | 2D spatial heatmap comparing FLS, MILS, and MFLS at t = 10 yr; asymmetric plume visible for moving models |

## Installation

The package is registered in the Julia General registry:

```julia
using Pkg
Pkg.add("GroundResponse")
```

Or in the Julia REPL package manager (`]`):

```
pkg> add GroundResponse
```

## Dependencies

### Library

| Package | Used for |
|---------|----------|
| [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl) | `expinti`, `erf`, Bessel functions in ILS / ICS / MILS / MFLS kernels |
| [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl) | Numerical integration in FLS, ICS, MFLS |
| [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) | Matrix operations (stdlib) |

### Scripts only

| Package | Used in |
|---------|---------|
| [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl) | `script_ground_models.jl`, `script_groundwater_advection.jl` |

## References

### Ground Thermal Response Models

- Ingersol, L. R. (1948). Theory of the ground pipe heat source for the heat pump. *ASHVE Journal Section, Heating, Piping and Air Conditioning*.
- Carslaw, H. S., & Jaeger, J. C. (1959). *Conduction of Heat in Solids* (2nd ed.). Oxford: Clarendon Press.
- Claesson, J., & Javed, S. (2011). An analytical method to calculate borehole fluid temperatures for time-scales from minutes to decades. *ASHRAE Transactions*, 117(PART 2), 279–288.
- Cimmino, M., & Bernier, M. (2014). A semi-analytical method to generate g-functions for geothermal bore fields. International Journal of Heat and Mass Transfer, 70, 641–650. https://doi.org/10.1016/j.ijheatmasstransfer.2013.11.037
- Guo, Y., Hu, X., Banks, J., & Liu, W. V. (2020). Considering buried depth in the moving finite line source model for vertical borehole heat exchangers — A new solution. *Energy and Buildings*, 214, 109859. https://doi.org/10.1016/j.enbuild.2020.109859
- Pasquier, P., & Lamarche, L. (2022). Analytic expressions for the moving infinite line source model. *Geothermics*, 103, 102413. https://doi.org/10.1016/j.geothermics.2022.102413

### Response Arrays

- Cimmino, M. (2015). The effects of borehole thermal resistances and fluid flow rate on the g-functions of geothermal bore fields. *International Journal of Heat and Mass Transfer*, 91, 1119–1127. https://doi.org/10.1016/j.ijheatmasstransfer.2015.08.041
- Dusseault, B., Pasquier, P., & Marcotte, D. (2018). A block matrix formulation for efficient g-function construction. *Renewable Energy*, 121, 249–260. https://doi.org/10.1016/j.renene.2017.12.092
