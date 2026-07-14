module GroundResponse

using SpecialFunctions
using QuadGK
using LinearAlgebra
using DSP
using PCHIPInterpolation: Interpolator

# Structure for ground response models
"""
    AbstractGroundModel

Abstract type for ground thermal response models. Subtype this to define a ground model that can be
passed to `ground_response`, `successive_flux`, and `bloc_matrix`. Only specify the parameters of 
the model; the g-function backends are defined separately and will be called by `ground_response`.
"""
abstract type AbstractGroundModel end

"""
Infinite line source (Ingersol, 1948).
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
Finite line source (Claesson & Javed, 2011; Cimmino & Bernier, 2014)
Parameters `H` [m], `D` [m], `ks` [W/mK], `Cs` [J/m³K], and the optional `nseg` [-].

`nseg` is the number of segments each borehole is divided into and selects the boundary condition of
the spatial superposition:
  - `nseg = 1` (default): uniform heat flux along each borehole — the conventional FLS (BC-I/BC-II).
  - `nseg > 1`: heat flux varies along the borehole to keep the whole field at a uniform wall
    temperature — the segment boundary condition BC-III (Cimmino & Bernier, 2014). The four-argument
    constructor `FLSModel(H, D, ks, Cs)` keeps `nseg = 1`, so existing code is unaffected.
"""
struct FLSModel <: AbstractGroundModel
    H::Float64
    D::Float64
    ks::Float64
    Cs::Float64
    nseg::Int
end
FLSModel(H::Real, D::Real, ks::Real, Cs::Real) = FLSModel(H, D, ks, Cs, 1)

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
include("short_term_ann.jl")                # short-term ANN transfer function — Pasquier et al. (2018)

include("spatial_superposition.jl")
include("borefield.jl")

# High-level interface
"""
    ground_response(t, rb, xy, m::AbstractGroundModel; bc = :II, solver = :successive, interp = true)

Compute the g-function for either one borehole or a borefield using the given ground model.
- **Single borehole** (`size(xy, 1) == 1`): evaluates the model at the borehole wall radius `rb`
  (for `bc = :II`; `bc = :I` routes through `uniform_flux`). Subject to `interp` like any other path.
- **Multiple boreholes**: applies spatial superposition selected by `bc` / `solver`.
  `t` must be a vector of time steps.

# Arguments
    - `t`: Time value or vector [s]
    - `rb`: Borehole radius [m]
    - `xy`: Borehole coordinates (nb × 2) [m] — can be built with `borefield(:rectangle, ...)`
    - `m`: Ground model (e.g. `FLSModel(150, 4, 3.0, 2e6)`)
# Keywords
    - `bc`: Boundary condition of the spatial superposition:
        - `:I`   — equal, uniform heat flux on every borehole → [`uniform_flux`](@ref);
        - `:II`  — uniform flux along each borehole, equal mean wall temperature (default);
        - `:III` — flux varies within each borehole for a uniform wall temperature
          ([`FLSModel`](@ref) only; requires `nseg > 1`).
    - `solver`: backend for the chosen `bc`:
        - `bc = :II`  → `:successive` ([`successive_flux`](@ref), default) or `:block`
          ([`bloc_matrix`](@ref));
        - `bc = :III` → `:block` ([`segment_response`](@ref), default) or `:marching`
          ([`segment_response_marching`](@ref)).
    - `interp` (default `true`): compute on an internal constant-step sub-sampling grid and
      PCHIP-interpolate to `t`. One uniform keyword, but its *role* differs by backend:
        - **temporal solvers** (`:successive`, `:marching`): a **correctness** requirement — their
          convolution assumes a constant step, so this is the only way a non-uniform `t` yields a
          valid result. It also bounds the cost to the ~one hundred sub-sample nodes, independent of
          `length(t)`.
        - **instantaneous / direct backends** (`:block`, `bc = :I`, single borehole): a
          **performance** approximation only — these are exact at any `t`, so `interp` merely trades a
          small interpolation error for far fewer evaluations on a large `t`.
      Pass `interp = false` to compute directly at the requested `t` (the temporal solvers then
      require a uniformly spaced `t`).
# Output
    - `g`: g-function for an impulse of 1 W/m [°Cm/W]
"""
function ground_response(t, rb::Real, xy::AbstractMatrix{<:Real}, m::AbstractGroundModel;
    bc::Symbol = :II, solver::Symbol = :successive, interp::Bool = true)
    bc === :I && return _apply_interp(t, interp, tt -> uniform_flux(tt, rb, xy, m))
    if bc === :II
        size(xy, 1) == 1 && return _apply_interp(t, interp, tt -> _borehole_response(tt, rb, m))
        solver === :block      && return _apply_interp(t, interp, tt -> bloc_matrix(tt, rb, xy, m))
        solver === :successive && return successive_flux(t, rb, xy, m; interp)
        throw(ArgumentError("bc = :II solver must be :successive or :block, got :$solver"))
    end
    throw(ArgumentError("bc = :$bc is not supported for $(nameof(typeof(m))); use :I or :II"))
end

# FLSModel additionally supports BC-III (segment superposition). The default `bc` mirrors the
# historical behaviour: `nseg > 1` selects BC-III, `nseg == 1` stays on BC-II.
function ground_response(t, rb::Real, xy::AbstractMatrix{<:Real}, m::FLSModel;
    bc::Symbol = m.nseg > 1 ? :III : :II, solver::Symbol = :successive, interp::Bool = true)
    bc === :I && return _apply_interp(t, interp, tt -> uniform_flux(tt, rb, xy, m))
    if bc === :II
        size(xy, 1) == 1 && return _apply_interp(t, interp, tt -> _borehole_response(tt, rb, m))
        solver === :block      && return _apply_interp(t, interp, tt -> bloc_matrix(tt, rb, xy, m))
        solver === :successive && return successive_flux(t, rb, xy, m; interp)
        throw(ArgumentError("bc = :II solver must be :successive or :block, got :$solver"))
    elseif bc === :III
        solver === :marching && return segment_response_marching(t, rb, xy, m; interp)
        return _apply_interp(t, interp, tt -> segment_response(tt, rb, xy, m))  # block (default BC-III)
    end
    throw(ArgumentError("bc must be :I, :II, or :III, got :$bc"))
end

_borehole_response(t, rb, m::ILSModel)  = ils(t, rb, m.ks, m.Cs)
_borehole_response(t, rb, m::ICSModel)  = ics(t, rb, rb, m.ks, m.Cs)
_borehole_response(t, rb, m::FLSModel)  = fls(t, rb, m.H, m.D, m.ks, m.Cs)
# The borehole-wall self response is the circumferential mean at r = rb (θ irrelevant inside).
_borehole_response(t, ::Real, m::MILSModel) = mils(t, m.rb, 0.0, m.rb, m.ks, m.Cs, m.Cf, m.vD)
_borehole_response(t, ::Real, m::MFLSModel) = mfls(t, m.rb, 0.0, m.H, m.rb, m.D, m.ks, m.Cs, m.Cf,
    m.vD)

# -------------------------------------------------------------------------
# Exports
# -------------------------------------------------------------------------

# Abstract type — downstream packages dispatch on this
export AbstractGroundModel
export ILSModel, ICSModel, FLSModel, MILSModel, MFLSModel

# Raw backends — available for direct use
export ils, ics, fls, mils, mfls
export short_term_response

# High-level interface
export ground_response

# Spatial superposition
export bloc_matrix, successive_flux, uniform_flux, segment_response, segment_response_marching

# Borefield layouts
export borefield_geometry
export borefield, borefield_rectangle, borefield_line, borefield_circle
export borefield_L, borefield_U, borefield_open_rectangle

end
