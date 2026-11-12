# GroundResponse.jl

[![CI](https://github.com/GHE-jl/GroundResponse.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/GHE-jl/GroundResponse.jl/actions/workflows/CI.yml)
[![Docs: dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GHE-jl.github.io/GroundResponse.jl/dev)
[![Docs: stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://GHE-jl.github.io/GroundResponse.jl/stable)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPLv3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

A Julia package providing analytical ground thermal response functions (g-functions) for single
boreholes and borehole fields. GroundResponse.jl is designed as a computational backbone: the raw
model kernels (`ils`, `fls`, …) evaluate uninterpolated responses at arbitrary time steps and radii,
leaving system-level modelling to downstream packages such as
[GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl).

Time discretisation *is* handled here, though, because it is intrinsic to the spatial
superposition: under BC-II and BC-III the flux of each borehole changes in time, and its load
history is superposed on a **uniform time grid** (`t = Δt, 2Δt, …`), so an arbitrary (e.g.
log-spaced) time vector would otherwise give a wrong answer. The high-level interface therefore
solves on internal uniform blocks and interpolates back to the requested times — see the `interp`
keyword below.

## Models

Five analytical ground response models are implemented, each as a raw function and as an
`AbstractGroundModel` subtype for the high-level interface.

| Model | Function | Struct | Parameters |
|-------|----------|--------|------------|
| Infinite line source (Ingersol, 1948) | `ils` | `ILSModel` | `ks`, `Cs` |
| Infinite cylindrical source (Carslaw & Jaeger, 1959) | `ics` | `ICSModel` | `rc`, `ks`, `Cs` |
| Finite line source (Claesson & Javed, 2011) | `fls` | `FLSModel` | `H`, `D`, `ks`, `Cs` |
| Moving infinite line source (Pasquier & Lamarche, 2022) | `mils` | `MILSModel` | `rb`, `ks`, `Cs`, `Cf`, `vD` |
| Moving finite line source (Guo et al., 2020) | `mfls` | `MFLSModel` | `H`, `rb`, `D`, `ks`, `Cs`, `Cf`, `vD` |

**Parameter legend** — `H`: borehole depth [m], `D`: buried depth [m], `rb`: borehole radius [m],
`rc`: cylinder radius [m], `ks`: ground thermal conductivity [W/mK], `Cs`: ground volumetric heat
capacity [J/m³K], `Cf`: groundwater volumetric heat capacity [J/m³K], `vD`: Darcy velocity [m/s].

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

## High-Level Interface

`ground_response` is a single entry point that dispatches over model type and borefield size:

```julia
g = ground_response(t, rb, xy, m::AbstractGroundModel; bc = :II, interp = true)
```

- **Single borehole** (`size(xy, 1) == 1`): evaluates the model kernel at the borehole wall.
- **Multiple boreholes**: computes the borefield g-function via the selected spatial superposition.

Keywords:

- **`bc`** — boundary condition of the spatial superposition:
  - `:I` — equal, uniform heat flux on every borehole (`bc_I`);
  - `:II` — uniform flux per borehole, equal mean wall temperature (`bc_II`, default);
  - `:III` — flux varies within each borehole for a uniform wall temperature (`bc_III`; `FLSModel`
    with `nseg > 1`, where it is the default).
- **`interp`** (default `true`) — solve on internal uniform blocks (hourly to one day, daily to one
  month, monthly to one year, quarterly to ten years, yearly to fifty years) and PCHIP-interpolate
  to `t` in log time. **One keyword, two roles**:
  - for **BC-II and BC-III** it is a *correctness* requirement on any `t` that is not
    `Δt, 2Δt, …, nΔt` (their load history needs a uniform grid), and it also bounds the cost to
    about two hundred solved times, independent of `length(t)`;
  - for **BC-I and a single borehole** it is a *pure performance* approximation — those are exact
    at any `t`, so `interp` just trades a small interpolation error for far fewer evaluations on a
    large `t`.

  Pass `interp = false` to compute directly at the requested `t` (BC-II and BC-III then require
  `t = Δt, 2Δt, …, nΔt`, and reject any other `t`).
- **`solver`** — deprecated and ignored since 2.2.0 (BC-II and BC-III are always solved by time
  marching).

```julia
using GroundResponse

ks, Cs = 3.0, 2.0e6          # Ground thermal conductivity [W/mK], heat capacity [J/m³K]
rb     = 0.076               # Borehole radius [m]
t      = 3600.0 .* exp10.(range(0, log10(8760*25), length=200))  # 1 h → 25 yr

# Single borehole — FLS model
m = FLSModel(150.0, 4.0, ks, Cs)           # H=150 m, D=4 m
g = ground_response(t, rb, [0.0 0.0], m)   # g-function [°Cm/W]

# Borefield — same model, spatial superposition applied automatically
xy = borefield(:rectangle, 3, 3, 5.0)      # 3×3 grid, 5 m spacing [m]
g  = ground_response(t, rb, xy, m)

# Moving model (with groundwater advection)
Cf, vD = 4.2e6, 1e-6
m_mfls = MFLSModel(150.0, rb, 4.0, ks, Cs, Cf, vD)
g = ground_response(t, rb, xy, m_mfls)
```

## Spatial Superposition

One function per boundary condition computes the borefield g-function from the single-borehole
responses. Each accepts either a precomputed `nt × nb × nb` g-matrix (low-level kernel) or borefield
parameters with any `AbstractGroundModel` (high-level, dispatching through `_borehole_response`):

| Function | BC | Solution | `t`-spacing |
|----------|----|----------|-------------|
| `bc_I` | I — same flux on every borehole (Guo et al., 2021) | sum, no solve | independent |
| `bc_II` | II — uniform flux per borehole, equal mean wall temperature | direct, time marching | `Δt, 2Δt, …` |
| `bc_III` | III — segments, uniform wall temperature (Cimmino & Bernier, 2014) | direct, time marching | `Δt, 2Δt, …` |

Under BC-II and BC-III the flux of each borehole (or segment) changes in time, so each time step
solves a small system for the flux increments and the common wall temperature, with the load
history of the previous steps on the right-hand side. This is the block-matrix system of Dusseault
et al. (2018) solved by forward substitution: same solution, no iteration, a fraction of the cost.
`bc_III` with `nseg = 1` is `bc_II`.

```julia
# High-level: pass the model directly (interp = true by default)
g = bc_I(t, rb, xy, MILSModel(rb, ks, Cs, Cf, vD))
g = bc_II(t, rb, xy, FLSModel(150.0, 4.0, ks, Cs))
g = bc_III(t, rb, xy, FLSModel(150.0, 4.0, ks, Cs, 8))

# Low-level: pass a precomputed 3D g-matrix on t = Δt, 2Δt, …
tu  = 2_592_000.0 .* (1:240)                      # monthly steps over 20 years [s]
r3D = borefield_geometry(xy, rb)[1]               # nb×nb distance matrix
g3D = fls(tu, r3D, 150.0, 4.0, ks, Cs)           # nt × nb × nb
g   = bc_II(g3D)
```

The solvers used before 2.2.0 are kept for comparison as `successive_flux_legacy` (iterative,
Nguyen & Pasquier, 2021), `bloc_matrix_legacy` and `segment_response_legacy` (block matrix, corrected
in 2.2.0 to include the load history) and `segment_response_marching_legacy`. Their former names
(`successive_flux`, `bloc_matrix`, `segment_response`, `segment_response_marching`, and
`uniform_flux` for `bc_I`) still work as deprecated aliases.

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

Subtype `AbstractGroundModel` and add a single `_borehole_response` method returning the pairwise
`nt × nb × nb` response array. `bc_I`, `bc_II` and the `ground_response` interface — including
`interp` sub-sampling — then work automatically:

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
| `script_spatial_superposition.jl` | BC-I / BC-II / BC-III on a 3×3 field; BC-III convergence with segments; g-function scaling from 1×1 to 5×5; computing time of `bc_II` / `bc_III` against the legacy solvers; `interp = true` against an hourly direct solve |
| `script_borefield_layouts.jl` | Subplot figure of all six layout configurations (~50 boreholes each); g-functions of two fields rotated from 0 to 90° against the groundwater flow (FLS and MFLS) |
| `script_ground_response.jl` | `ground_response` dispatch across all models; `@btime` benchmark for single borehole and borefield sizes |
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
| [DSP.jl](https://github.com/JuliaDSP/DSP.jl) | `conv` for the legacy successive flux solver (`successive_flux_legacy`) |
| [PCHIPInterpolation.jl](https://github.com/gerlero/PCHIPInterpolation.jl) | Interpolating the sub-sampled g-function back to the requested times (`interp`) |
| [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) | Matrix operations (stdlib) |

### Scripts only

| Package | Used in |
|---------|---------|
| [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl) | All visualisation scripts |
| [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl) | `script_ground_response.jl` |

## References

### Ground Thermal Response Models

- Ingersol, L. R. (1948). Theory of the ground pipe heat source for the heat pump. *ASHVE Journal Section, Heating, Piping and Air Conditioning*.
- Carslaw, H. S., & Jaeger, J. C. (1959). *Conduction of Heat in Solids* (2nd ed.). Oxford: Clarendon Press.
- Claesson, J., & Javed, S. (2011). An analytical method to calculate borehole fluid temperatures for time-scales from minutes to decades. *ASHRAE Transactions*, 117(PART 2), 279–288.
- Cimmino, M., & Bernier, M. (2014). A semi-analytical method to generate g-functions for geothermal bore fields. International Journal of Heat and Mass Transfer, 70, 641–650. https://doi.org/10.1016/j.ijheatmasstransfer.2013.11.037
- Guo, Y., Hu, X., Banks, J., & Liu, W. V. (2020). Considering buried depth in the moving finite line source model for vertical borehole heat exchangers — A new solution. *Energy and Buildings*, 214, 109859. https://doi.org/10.1016/j.enbuild.2020.109859
- Guo, Y., Hu, X., Banks, J., & Liu, W. V. (2021). Considering buried depth for vertical borehole heat exchangers in a borehole field with groundwater flow — An extended solution. Energy and Buildings, 235, 110722. https://doi.org/10.1016/j.enbuild.2021.110722
- Pasquier, P., & Lamarche, L. (2022). Analytic expressions for the moving infinite line source model. *Geothermics*, 103, 102413. https://doi.org/10.1016/j.geothermics.2022.102413

### Spatial Superposition

- Cimmino, M., Bernier, M., & Adams, F. (2013). A contribution towards the determination of g-functions using the finite line source. *Applied Thermal Engineering*, 51(1–2), 401–412. https://doi.org/10.1016/j.applthermaleng.2012.07.044
- Cimmino, M. (2018). Fast calculation of the g-functions of geothermal borehole fields using similarities in the evaluation of the finite line source solution. *Journal of Building Performance Simulation*, 11(6), 655–668. https://doi.org/10.1080/19401493.2017.1423390
- Dusseault, B., Pasquier, P., & Marcotte, D. (2018). A block matrix formulation for efficient g-function construction. *Renewable Energy*, 121, 249–260. https://doi.org/10.1016/j.renene.2017.12.092
- Nguyen, A., & Pasquier, P. (2021). A successive flux estimation method for rapid g-function construction of small to large-scale ground heat exchanger. *Renewable Energy*, 165, 359–368. https://doi.org/10.1016/j.renene.2020.10.074
