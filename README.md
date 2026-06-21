# GroundResponse.jl

A Julia package providing analytical ground response functions (g-functions) for single boreholes and borehole fields. GroundResponse.jl serves as a backbone to [GroundHeatExchanger.jl](https://github.com/gabriel-dion/GroundHeatExchanger.jl), computing uninterpolated thermal responses at arbitrary time steps and radii and leaving temporal sampling and interpolation to downstream packages.

## Models

Five analytical ground response models are implemented both as raw functions and as `AbstractGroundModel` subtypes for the high-level interface.

⚠️A short-term ANN model is also included, but not yet fully implemented.

| Model | Function | Struct | Parameters |
|-------|----------|--------|------------|
| Infinite line source (Ingersol, 1948) | `ils` | `ILSModel` | `ks`, `Cs` |
| Infinite cylindrical source (Carslaw & Jaeger, 1959) | `ics` | `ICSModel` | `rc`, `ks`, `Cs` |
| Finite line source (Claesson & Javed, 2011) | `fls` | `FLSModel` | `H`, `D`, `ks`, `Cs` |
| Moving infinite line source (Pasquier & Lamarche, 2022) | `mils` | `MILSModel` | `rb`, `ks`, `Cs`, `Cf`, `vD` |
| Moving finite line source (Guo et al., 2020) | `mfls` | `MFLSModel` | `H`, `rb`, `D`, `ks`, `Cs`, `Cf`, `vD` |
| Short-term ANN (Pasquier et al., 2018) | `gST_ANN` | — | — |

**Parameters**: `H` — borehole depth [m], `D` — buried depth [m], `rb` — borehole radius [m], `rc` — cylinder radius [m], `ks` — ground thermal conductivity [W/mK], `Cs` — ground volumetric heat capacity [J/m³K], `Cf` — groundwater volumetric heat capacity [J/m³K], `vD` — Darcy velocity [m/s]. The moving models require `vD > 0`; use a small value (e.g. `1e-12`) for near-impervious conditions or the corresponding non-moving model instead.

The MILS and MFLS models are directional: both implement a single unified kernel that switches between the circumferential-average form (inside/at the borehole, `r ≤ rb`) and the direction-dependent form (outside, `r > rb`). Groundwater flow is assumed along the positive x-axis.

## Dispatch

All models are available as direct functions with multiple dispatch over time `t` and the spatial argument. The output shape is determined entirely by the input types, with no runtime branching on values.

### ILS, ICS, FLS — scalar radius `r`

| `t` | `r` | Output | Notes |
|-----|-----|--------|-------|
| `Real` | `Real` | scalar | Single time step, single radius |
| `AbstractVector` | `Real` | `nt`-vector | Time series at one radius |
| `Real` | `AbstractVector` | `nr`-vector | Radial profile at one time |
| `AbstractVector` | `AbstractVector` | `nt × nr` matrix | Time series at multiple radii |
| `Real` | `AbstractMatrix` (nb×nb) | `nb × nb` matrix | Single time step, borefield radius matrix |
| `AbstractVector` | `AbstractMatrix` (nb×nb) | `nt × nb × nb` 3D array | Full borefield g-matrix |

The `nb×nb` pairwise radius matrix for a borefield is computed by `borefield_radius(xy, rb)`, which also returns unique radii, indices, and angle vectors.

### MILS, MFLS — directional coordinates `xy`

| `t` | `xy` | Output | Notes |
|-----|------|--------|-------|
| `Real` | `AbstractVector` ([x, y]) | scalar | Single time step, single point |
| `AbstractVector` | `AbstractVector` ([x, y]) | `nt`-vector | Time series at one point |
| `Real` | `AbstractMatrix` (nb×2) | `nb × nb` matrix | Single time step, full borefield |
| `AbstractVector` | `AbstractMatrix` (nb×2) | `nt × nb × nb` 3D array | Full borefield g-matrix |

For the matrix overloads, `xy` is the `nb×2` matrix of absolute borehole coordinates (one borehole per row). The function internally computes all pairwise displacement vectors from each source borehole to each target; the diagonal (self-response) uses displacement `[0, 0]`, which triggers the inside-borehole branch of the kernel. The `nb×2` coordinate matrix can be built with `borefield_xy(...)`.

## High-Level Interface

`ground_response` provides a single entry point that dispatches over the model type and borefield size:

```julia
g = ground_response(t, rb, xy, m::AbstractGroundModel)
```

- **Single borehole** (`size(xy, 1) == 1`): calls the model directly at the borehole wall.
- **Multiple boreholes**: computes the spatial superposition via `successive_flux`.

```julia
using GroundResponse

ks, Cs = 3.0, 2e6             # Ground properties [W/mK], [J/m³K]
rb     = 0.076                # Borehole radius [m]
t      = 3600.0 .* (1:8760)   # Hourly steps over one year [s]

# Single borehole
m = FLSModel(150.0, 4.0, ks, Cs)
g = ground_response(t, rb, [0.0 0.0], m)   # g-function [°Cm/W]

# Borefield
xy = borefield_xy(3, 3, 5.0)               # 3×3 grid, 5 m spacing [m]
g  = ground_response(t, rb, xy, m)

# Moving model (with groundwater)
Cf, vD = 4.2e6, 1e-6
m_mfls = MFLSModel(150.0, rb, 4.0, ks, Cs, Cf, vD)
g = ground_response(t, rb, xy, m_mfls)
```

## Spatial Superposition

Two methods aggregate single-borehole g-functions into a borefield g-function. Both accept either a precomputed 3D g-matrix `(nt × nb × nb)` or borefield parameters with any `AbstractGroundModel`:

- `successive_flux` — iterative successive flux method (Nguyen & Pasquier, 2021). Converges rapidly for small to large fields.
- `bloc_matrix` — block matrix formulation (Dusseault et al., 2018). Direct solve, suitable when the full system matrix is needed.

```julia
g = successive_flux(t, rb, xy, FLSModel(150.0, 4.0, ks, Cs))
g = bloc_matrix(t, rb, xy, MILSModel(rb, ks, Cs, Cf, vD))

# Or pass a precomputed 3D g-matrix directly
g3D = fls(t, borefield_radius(xy, rb)[1], 150.0, 4.0, ks, Cs)  # nt × nb × nb
g   = successive_flux(g3D)
```

## Extending with Custom Models

Subtype `AbstractGroundModel` and dispatch `successive_flux` and `bloc_matrix` on the new type:

```julia
struct MyModel <: AbstractGroundModel
    ks::Float64
    Cs::Float64
end

function GroundResponse.successive_flux(t, rb, xy, m::MyModel)
    r = borefield_radius(xy, rb)[1]          # nb×nb radius matrix
    return successive_flux(my_gfunc(t, r, m.ks, m.Cs))
end
```

## Utility Functions

- `borefield_xy(nx, ny, D)` — `nb×2` coordinate matrix for a rectangular `nx×ny` grid with spacing `D` [m].
- `borefield_radius(xy, rb)` — `nb×nb` pairwise radius matrix, unique radii, indices, and angle vector for a borefield.

## Installation

The package is not yet registered. Install directly from the repository:

```julia
using Pkg
Pkg.add(url="https://github.com/gabriel-dion/GroundResponse.jl")
```

Or in the Julia REPL package manager (`]`):

```
pkg> add https://github.com/gabriel-dion/GroundResponse.jl
```

## Dependencies

- [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl)
- [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl)
- [DSP.jl](https://github.com/JuliaDSP/DSP.jl)
- [PCHIPInterpolation.jl](https://github.com/gerlero/PCHIPInterpolation.jl) (used by `gST_ANN`)
- [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) (standard library)

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

## License

See [LICENSE.txt](LICENSE.txt).