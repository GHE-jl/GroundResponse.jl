module GroundResponse

using SpecialFunctions
using QuadGK
using LinearAlgebra
using DSP

# Structure for ground response models
"""
    AbstractGroundModel

Abstract type for ground thermal response models. Subtype this to define a ground model that can be
passed to `ground_response`, `successive_flux`, and `bloc_matrix`. Only specify the parameters of 
the model; the g-function backends are defined separately and will be called by `ground_response`.
"""
abstract type AbstractGroundModel end

"""
Infinite line source (Ingersol, 1954)
Parameters `ks` [W/mK], `Cs` [J/m³K].
"""
struct ILSModel <: AbstractGroundModel
    ks::Float64
    Cs::Float64
end

"""
Infinite cylindrical source (Carslaw & Jaeger, 1959):
Parameters `rc` [m], `ks` [W/mK], `Cs` [J/m³K].
"""
struct ICSModel <: AbstractGroundModel
    rc::Float64
    ks::Float64
    Cs::Float64
end

"""
Finite line source (Claesson & Javed, 2011)
Parameters `H` [m], `D` [m], `ks` [W/mK], `Cs` [J/m³K].
"""
struct FLSModel <: AbstractGroundModel
    H::Float64
    D::Float64
    ks::Float64
    Cs::Float64
end

"""
Moving infinite line source (Pasquier & Lamarche, 2022)
Parameters `rb` [m], `ks` [W/mK], `Cs` [J/m³K], `Cf` [J/m³K], `vD` [m/s].
"""
struct MILSModel <: AbstractGroundModel
    rb::Float64
    ks::Float64
    Cs::Float64
    Cf::Float64
    vD::Float64
end

"""
Moving finite line source (Guo et al., 2020)
Parameters `H` [m], `rb` [m], `D` [m], `ks` [W/mK], `Cs` [J/m³K], `Cf` [J/m³K], `vD` [m/s].
"""
struct MFLSModel <: AbstractGroundModel
    H::Float64
    rb::Float64
    D::Float64
    ks::Float64
    Cs::Float64
    Cf::Float64
    vD::Float64
end

include("infinite_line_source.jl")          # ils  — Ingersol (1948)
include("infinite_cylindrical_source.jl")   # ics  — Carslaw & Jaeger (1959)
include("finite_line_source.jl")            # fls  — Claesson & Javed (2011)
include("moving_infinite_line_source.jl")   # mils — Pasquier & Lamarche (2022)
include("moving_finite_line_source.jl")     # mfls — Guo et al. (2020)
include("g_short_term.jl")                  # gST_ANN — Pasquier et al. (2018)

include("spatial_superposition.jl")
include("utils.jl")

# High-level interface
"""
    ground_response(t, rb, xy, m::AbstractGroundModel)

Compute the g-function for either one borehole or a borefield using the given ground model.
- **Single borehole** (`size(xy, 1) == 1`): evaluates the model directly at the
  borehole wall radius `rb`.
- **Multiple boreholes**: applies spatial superposition via `successive_flux`.
  `t` must be a vector of time steps (or a scalar).
# Arguments
    - `t`: Time value or vector [s]
    - `rb`: Borehole radius [m]
    - `xy`: Borehole coordinates (nb × 2) [m] — can be built with `borefield_xy`
    - `m`: Ground model (e.g. `FLSModel(150, 4, 3.0, 2e6)`)
# Output
    - `g`: g-function for an impulse of 1 W/m [°Cm/W]
"""
function ground_response(t, rb::Real, xy::AbstractMatrix{<:Real}, m::AbstractGroundModel)
    if size(xy, 1) == 1
        return _borehole_response(t, rb, m)
    else
        return successive_flux(t, rb, xy, m)
    end
end

_borehole_response(t, rb, m::ILSModel)  = ils(t, rb, m.ks, m.Cs)
_borehole_response(t, rb, m::ICSModel)  = ics(t, rb, rb, m.ks, m.Cs)
_borehole_response(t, rb, m::FLSModel)  = fls(t, rb, m.H, m.D, m.ks, m.Cs)
_borehole_response(t, ::Real, m::MILSModel) = mils(t, [0.0, 0.0], m.rb, m.ks, m.Cs, m.Cf, m.vD)
_borehole_response(t, ::Real, m::MFLSModel) = mfls(t, [0.0, 0.0], m.H, m.rb, m.D, m.ks, m.Cs, m.Cf, m.vD)

# -------------------------------------------------------------------------
# Exports
# -------------------------------------------------------------------------

# Abstract type — downstream packages dispatch on this
export AbstractGroundModel
export ILSModel, ICSModel, FLSModel, MILSModel, MFLSModel

# Raw backends — available for direct use
export ils, ics, fls, mils, mfls

# High-level interface
export ground_response

# Spatial superposition
export bloc_matrix, successive_flux

# Utilities
export borefield_radius, borefield_xy

end
