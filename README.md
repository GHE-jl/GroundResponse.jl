# GroundResponse.jl

A Julia package providing analytical ground thermal response functions (g-functions) for single
boreholes and borehole fields. GroundResponse.jl is designed as a computational backbone: it
evaluates uninterpolated thermal responses at arbitrary time steps and radii, leaving temporal
sampling, interpolation, and system-level modelling to downstream packages such as
[GroundHeatExchanger.jl](https://github.com/GeothermalJL/GroundHeatExchanger.jl).

---

## Models

Five analytical ground response models are implemented, each as a raw function and as an
`AbstractGroundModel` subtype for the high-level interface. A short-term ANN model is also
included but is not yet fully integrated (⚠️ see `gST_ANN`).

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

---

## Multiple Dispatch

All functions are overloaded on the type of `t` and the spatial argument, with no runtime
branching. The output shape is determined entirely by the input types.

### ILS, ICS, FLS — scalar radius `r`

| `t` | `r` | Output | Notes |
|-----|-----|--------|-------|
| `Real` | `Real` | scalar | Single time step, single radius |
| `AbstractVector` | `Real` | `nt`-vector | Time series at one radius |
| `Real` | `AbstractVector` | `nr`-vector | Radial profile at one time step |
| `AbstractVector` | `AbstractVector` | `nt × nr` matrix | Time series at multiple radii |
| `Real` | `AbstractMatrix` (nb×nb) | `nb × nb` matrix | Single time step, borefield radius matrix |
| `AbstractVector` | `AbstractMatrix` (nb×nb) | `nt × nb × nb` 3D array | Full borefield g-matrix |

The `nb×nb` pairwise radius matrix is produced by `borefield_radius(xy, rb)`, which also returns
the unique radii, index map, and azimuth angles.

### MILS, MFLS — directional coordinates `xy`

| `t` | `xy` | Output | Notes |
|-----|------|--------|-------|
| `Real` | `AbstractVector` ([x, y]) | scalar | Single time step, single point |
| `AbstractVector` | `AbstractVector` ([x, y]) | `nt`-vector | Time series at one point |
| `Real` | `AbstractMatrix` (nb×2) | `nb × nb` matrix | Single time step, full borefield |
| `AbstractVector` | `AbstractMatrix` (nb×2) | `nt × nb × nb` 3D array | Full borefield g-matrix |

For the matrix overloads, `xy` is the `nb×2` matrix of borehole coordinates. The kernel
internally computes pairwise displacements; the diagonal (self-response) uses `[0, 0]`, which
triggers the inside-borehole branch.

---

## High-Level Interface

`ground_response` is a single entry point that dispatches over model type and borefield size:

```julia
g = ground_response(t, rb, xy, m::AbstractGroundModel)
```

- **Single borehole** (`size(xy, 1) == 1`): calls the model kernel directly at the borehole wall.
- **Multiple boreholes**: computes the borefield g-function via `successive_flux`.

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

---

## Spatial Superposition

Two methods compute the borefield g-function from single-borehole responses. Both accept either a
precomputed `nt × nb × nb` g-matrix or borefield parameters with any `AbstractGroundModel`:

- **`successive_flux`** — iterative successive flux method (Nguyen & Pasquier, 2021). Fast and
  memory-efficient for small to large fields.
- **`bloc_matrix`** — block matrix formulation (Dusseault et al., 2018). Direct linear solve,
  more expensive but unconditionally stable.

```julia
# High-level: pass model directly
g = successive_flux(t, rb, xy, FLSModel(150.0, 4.0, ks, Cs))
g = bloc_matrix(t, rb, xy, MILSModel(rb, ks, Cs, Cf, vD))

# Low-level: pass precomputed 3D g-matrix
r3D = borefield_radius(xy, rb)[1]                 # nb×nb radius matrix
g3D = fls(t, r3D, 150.0, 4.0, ks, Cs)            # nt × nb × nb
g   = successive_flux(g3D)
```

---

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

All functions return an `nb×2` matrix of borehole coordinates `[x y]`.

`borefield_radius(xy, rb)` returns `(r, rᵥ, rᵤ, rᵢ, θ, nb)`: the pairwise radius matrix,
its flattened vector, unique values, index map, azimuth angles, and borehole count.

---

## Extending with Custom Models

Subtype `AbstractGroundModel` and add `successive_flux` / `bloc_matrix` overloads:

```julia
struct MyModel <: AbstractGroundModel
    ks::Float64
    Cs::Float64
end

function GroundResponse.successive_flux(t, rb, xy, m::MyModel)
    r, = borefield_radius(xy, rb)
    return successive_flux(my_gfunc(t, r, m.ks, m.Cs))
end

function GroundResponse.bloc_matrix(t, rb, xy, m::MyModel)
    r, = borefield_radius(xy, rb)
    return bloc_matrix(my_gfunc(t, r, m.ks, m.Cs))
end
```

The new model is then usable with `ground_response` without any changes to the core package.

---

## Scripts

Example scripts are in [`script/`](script/). Run them from the package root after activating the
project environment (`julia --project`):

| Script | Description |
|--------|-------------|
| [`script_ground_models.jl`](script/script_ground_models.jl) | Single-borehole g-function comparison: ILS, ICS, FLS, MILS, MFLS on a log-spaced time axis |
| [`script_spatial_superposition.jl`](script/script_spatial_superposition.jl) | Bloc matrix vs successive flux agreement (3×3 FLS); g-function scaling from 1×1 to 5×5 |
| [`script_borefield_layouts.jl`](script/script_borefield_layouts.jl) | Subplot figure of all six layout configurations (~50 boreholes each) |
| [`script_ground_response.jl`](script/script_ground_response.jl) | `ground_response` dispatch across all models; `@btime` benchmark for single borehole and borefield sizes |
| [`script_groundwater_advection.jl`](script/script_groundwater_advection.jl) | 2D spatial heatmap comparing FLS, MILS, and MFLS at t = 10 yr; asymmetric plume visible for moving models |
| [`script_short-term.jl`](script/script_short-term.jl) | Short-term ANN model (⚠️ not yet fully integrated) |

---

## Installation

The package is not yet registered. Install directly from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/GeothermalJL/GroundResponse.jl")
```

Or in the Julia REPL package manager (`]`):

```
pkg> add https://github.com/GeothermalJL/GroundResponse.jl
```

---

## Dependencies

### Library

| Package | Used for |
|---------|----------|
| [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl) | `expinti`, `erf`, Bessel functions in ILS / ICS / MILS / MFLS kernels |
| [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl) | Numerical integration in FLS, ICS, MFLS |
| [DSP.jl](https://github.com/JuliaDSP/DSP.jl) | `conv` for the successive flux iterative solver |
| [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) | Matrix operations (stdlib) |

### Scripts only

| Package | Used in |
|---------|---------|
| [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl) | All visualisation scripts |
| [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl) | `script_ground_response.jl` |

---

## References

### Ground Thermal Response Models

- Ingersol, L. R. (1948). Theory of the ground pipe heat source for the heat pump. *ASHVE Journal Section, Heating, Piping and Air Conditioning*.
- Carslaw, H. S., & Jaeger, J. C. (1959). *Conduction of Heat in Solids* (2nd ed.). Oxford: Clarendon Press.
- Claesson, J., & Javed, S. (2011). An analytical method to calculate borehole fluid temperatures for time-scales from minutes to decades. *ASHRAE Transactions*, 117(PART 2), 279–288.
- Guo, Y., Hu, X., Banks, J., & Liu, W. V. (2020). Considering buried depth in the moving finite line source model for vertical borehole heat exchangers — A new solution. *Energy and Buildings*, 214, 109859. https://doi.org/10.1016/j.enbuild.2020.109859
- Pasquier, P., & Lamarche, L. (2022). Analytic expressions for the moving infinite line source model. *Geothermics*, 103, 102413. https://doi.org/10.1016/j.geothermics.2022.102413
- Pasquier, P., Zarrella, A., & Labib, R. (2018). Application of artificial neural networks to near-instant construction of short-term g-functions. *Applied Thermal Engineering*. https://doi.org/10.1016/j.applthermaleng.2018.04.078

### Spatial Superposition

- Dusseault, B., Pasquier, P., & Marcotte, D. (2018). A block matrix formulation for efficient g-function construction. *Renewable Energy*, 121, 249–260. https://doi.org/10.1016/j.renene.2017.12.092
- Nguyen, A., & Pasquier, P. (2021). A successive flux estimation method for rapid g-function construction of small to large-scale ground heat exchanger. *Renewable Energy*, 165, 359–368. https://doi.org/10.1016/j.renene.2020.10.074

---

## License

See [LICENSE.txt](LICENSE.txt).
