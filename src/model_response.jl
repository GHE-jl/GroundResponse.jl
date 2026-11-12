# This file holds everything that connects a ground model to the rest of the package:
#   - `AbstractGroundModel`: type and the concrete model structs (their parameters);
#   - `borehole_response`: the single per-model extension point that turns a model into the
#       pairwise borehole response array the spatial-superposition backends consume;
#   - `ground_response`: the high-level dispatcher over the boundary conditions.

"""
    AbstractGroundModel

Abstract type for ground thermal response models. Subtype this to define a ground model that can be
passed to `ground_response`, and the spatial superposition backends. Only specify the parameters of
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
Parameters `H` [m], `D` [m], `ks` [W/mK], `Cs` [J/m³K], and optional `nseg` [-].

`nseg` is the number of segments each borehole is divided into and selects the boundary condition of
the spatial superposition:
    - `nseg = 1` (default): uniform heat flux along each borehole, the conventional FLS (BC-I/BC-II)
    - `nseg > 1`: heat flux varies along the borehole to keep the whole field at a uniform wall
        temperature — the segment boundary condition BC-III (Cimmino & Bernier, 2014)
"""
struct FLSModel <: AbstractGroundModel
    H::Float64
    D::Float64
    ks::Float64
    Cs::Float64
    nseg::Int
end
FLSModel(H::Real, D::Real, ks::Real, Cs::Real) = FLSModel(H, D, ks, Cs, 1) # Default nseg = 1

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

"""
    _check_time(t)

Reject non-positive requested times. Every kernel carries the Fourier number `αt/r²` in a
denominator, so `t ≤ 0` gives a singular zero-time response (`Inf`/`NaN`, or (worse) a silently
clamped value under `interp = true`).
# Arguments
    - `t`: Time value or vector [s]
"""
_time_error() = throw(ArgumentError(
    "all requested times must be strictly positive (t > 0); a zero or negative time gives a " *
    "singular zero-time response"))
_check_time(t::Real) = t > 0 ? nothing : _time_error()
_check_time(t::AbstractArray{<:Real}) = all(>(0), t) ? nothing : _time_error()

"""
    _borehole_response(t, rb, xy, m::AbstractGroundModel)

Borefield response under any ground model — the single per-model extension point. To add a new
model, subtype [`AbstractGroundModel`](@ref) and add one `borehole_response` method.
- **Multiple boreholes**: the pairwise response array `g[k, i, j]` — the response at borehole `i` to
    a unit flux at borehole `j` at time `t[k]` — sized `nt × nb × nb` (2D `nb × nb` for a scalar
    `t`).
- **Single borehole** (`size(xy, 1) == 1`): the borehole self response at the wall (`r = rb`),
    independent of the borehole's position, returned as a vector over `t` (scalar for a scalar `t`).

The isotropic models (ILS/ICS/FLS) key their geometry on the passed borehole radius `rb`; the moving
models (MILS/MFLS) carry their own `m.rb` and additionally need the flow-relative angle matrix
(their response is asymmetric under groundwater flow).
# Arguments
    - `t`: Time value or vector [s]
    - `rb`: Borehole radius [m]
    - `xy`: Borehole coordinates (nb × 2) [m] — can be built with the `borefield` function
    - `m`: Ground model (e.g. `FLSModel(150, 4, 3.0, 2e6)`)
"""
function _borehole_response(t, rb, xy, m::ILSModel)
    r = size(xy, 1) == 1 ? rb : borefield_geometry(xy, rb)[1]
    return ils(t, r, m.ks, m.Cs)
end
function _borehole_response(t, rb, xy, m::ICSModel)
    r = size(xy, 1) == 1 ? rb : borefield_geometry(xy, rb)[1]
    return ics(t, r, rb, m.ks, m.Cs)
end
function _borehole_response(t, rb, xy, m::FLSModel)
    r = size(xy, 1) == 1 ? rb : borefield_geometry(xy, rb)[1]
    return fls(t, r, m.H, m.D, m.ks, m.Cs)
end
function _borehole_response(t, rb, xy, m::MILSModel)
    rb == m.rb || throw(ArgumentError(
        "borehole radius rb ($rb) must match the model's rb ($(m.rb))"))
    size(xy, 1) == 1 && return mils(t, m.rb, 0.0, m.rb, m.ks, m.Cs, m.Cf, m.vD)
    r, θ = borefield_geometry(xy, m.rb)
    return mils(t, r, θ, m.rb, m.ks, m.Cs, m.Cf, m.vD)
end
function _borehole_response(t, rb, xy, m::MFLSModel)
    rb == m.rb || throw(ArgumentError(
        "borehole radius rb ($rb) must match the model's rb ($(m.rb))"))
    size(xy, 1) == 1 && return mfls(t, m.rb, 0.0, m.H, m.rb, m.D, m.ks, m.Cs, m.Cf, m.vD)
    r, θ = borefield_geometry(xy, m.rb)
    return mfls(t, r, θ, m.H, m.rb, m.D, m.ks, m.Cs, m.Cf, m.vD)
end

"""
    ground_response(t, rb, xy, m::AbstractGroundModel; bc = :II, interp = true)
    ground_response(t, rb, xy, m::FLSModel; bc = m.nseg > 1 ? :III : :II, interp = true)

Compute the g-function of one borehole or of a borefield with the given ground model.
- **Single borehole** (`size(xy, 1) == 1`): the model evaluated at the borehole wall radius `rb`.
- **Multiple boreholes**: spatial superposition under the boundary condition `bc`, by
    [`bc_I`](@ref), [`bc_II`](@ref) or [`bc_III`](@ref).
# Arguments
    - `t`: Time value or vector [s]
    - `rb`: Borehole radius [m]
    - `xy`: Borehole coordinates (nb × 2) [m], e.g. from [`borefield`](@ref)
    - `m`: Ground model (e.g. `FLSModel(150, 4, 3.0, 2e6)`)
# Keywords
    - `bc`: Boundary condition of the spatial superposition:
        - `:I` - equal, uniform heat flux on every borehole ([`bc_I`](@ref));
        - `:II` - uniform flux along each borehole, equal mean wall temperature ([`bc_II`](@ref),
            default);
        - `:III` - flux free to vary between the segments of each borehole, uniform wall
            temperature ([`bc_III`](@ref), [`FLSModel`](@ref) only). The default for an
            [`FLSModel`](@ref) with `nseg > 1`.
    - `interp` (default `true`): compute on the internal sub-sampling grid (uniform blocks from one
        hour to decades) and PCHIP-interpolate to `t`. Required by BC-II and BC-III for any `t` that
        is not `Δt, 2Δt, …, nΔt`. For BC-I and a single borehole it is a performance choice only.
        Pass `interp = false` to compute directly at `t`.
# Output
    - `g`: g-function for 1 W/m on every borehole [°C·m/W]
"""
function ground_response(t, rb::Real, xy::AbstractMatrix{<:Real}, m::AbstractGroundModel;
    bc::Symbol = :II, interp::Bool = true)
    _check_time(t)

    # Single borehole: the borefield response reduces to the borehole self response
    size(xy, 1) == 1 && return _apply_interp(t, interp, tt -> _borehole_response(tt, rb, xy, m))

    bc === :I  && return bc_I(t, rb, xy, m; interp)
    bc === :II && return bc_II(t, rb, xy, m; interp)

    # BC-III (segment superposition) is only available for an FLSModel; see the method below
    throw(ArgumentError("bc = :$bc is not supported for $(nameof(typeof(m))); use :I or :II"))
end
function ground_response(t, rb::Real, xy::AbstractMatrix{<:Real}, m::FLSModel;
    bc::Symbol = m.nseg > 1 ? :III : :II, interp::Bool = true)
    _check_time(t)

    # Single borehole. Under a uniform flux (BC-I/BC-II) the borefield response is just the
    # borehole self response. Under BC-III the flux varies along the borehole, so it goes to the
    # segment solver.
    size(xy, 1) == 1 && bc !== :III &&
        return _apply_interp(t, interp, tt -> _borehole_response(tt, rb, xy, m))

    bc === :I   && return bc_I(t, rb, xy, m; interp)
    bc === :II  && return bc_II(t, rb, xy, m; interp)
    bc === :III && return bc_III(t, rb, xy, m; interp)

    throw(ArgumentError("bc = :$bc is not supported; use :I, :II, or :III"))
end
