# This file holds everything that connects a ground model to the borehole responses:
#   - `AbstractGroundModel`: type and the concrete model structs (their parameters);
#   - `gmatrix`: the single per-model extension point that turns a model into the pairwise
#       borehole (or segment) response array. The spatial superposition of these arrays under a
#       boundary condition (BC-I to BC-IV) is done by GroundHeatExchanger.jl.

"""
    AbstractGroundModel

Abstract type for ground thermal response models. Subtype this to define a ground model and add a
[`gmatrix`](@ref) method for it. Only specify the parameters of the model; the spatial
superposition of the responses is done by GroundHeatExchanger.jl.
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
Infinite cylindrical source (Carslaw & Jaeger, 1959), with the borehole radius as the cylinder
radius. Parameters `ks` [W/mK], `Cs` [J/m³K].
"""
struct ICSModel <: AbstractGroundModel
    ks::Float64
    Cs::Float64
end

"""
Finite line source (Claesson & Javed, 2011; Cimmino & Bernier, 2014)
Parameters `H` [m], `D` [m], `ks` [W/mK], `Cs` [J/m³K].

The borehole can be split into segments with the second method of [`gmatrix`](@ref).
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

"""
    _check_time(t)

Reject non-positive requested times. Every kernel carries the Fourier number `αt/r²` in a
denominator, so `t ≤ 0` gives a singular zero-time response (`Inf`/`NaN`).
# Arguments
    - `t`: Time value or vector [s]
"""
_time_error() = throw(ArgumentError(
    "all requested times must be strictly positive (t > 0); a zero or negative time gives a " *
    "singular zero-time response"))
_check_time(t::Real) = t > 0 ? nothing : _time_error()
_check_time(t::AbstractArray{<:Real}) = all(>(0), t) ? nothing : _time_error()

"""
    gmatrix(t, rb, xy, m::AbstractGroundModel)
    gmatrix(t, rb, xy, m::FLSModel, nseg)

Response array of a borefield under the ground model `m`, the `G` of the block-matrix formulation
(Dusseault et al., 2018). To add a new model, subtype [`AbstractGroundModel`](@ref) and add one
`gmatrix` method.
- **Multiple boreholes**: `g[k, i, j]` is the response at borehole `i` to a unit flux (1 W/m) at
    borehole `j` at time `t[k]`, sized `nt × nb × nb` (`nb × nb` for a scalar `t`).
- **Single borehole** (`size(xy, 1) == 1`): the borehole self response at the wall (`r = rb`),
    independent of the borehole's position, returned as a vector over `t` (scalar for a scalar
    `t`).
- **Segments** (second method, finite line source only): each borehole is split into `nseg`
    segments of equal length, numbered from the top. `g[k, p, q]` is the mean response over segment
    `p` to a unit flux on segment `q`, sized `nt × nS × nS` with `nS = nb⋅nseg` and the segments
    ordered borehole-major (`p = (borehole − 1)⋅nseg + segment`). This is the array of the BC-III
    and BC-IV superpositions (Cimmino & Bernier, 2014; Cimmino, 2015).

The isotropic models (ILS/ICS/FLS) key their geometry on the passed borehole radius `rb`; the moving
models (MILS/MFLS) carry their own `m.rb` and additionally need the flow-relative angle matrix
(their response is asymmetric under groundwater flow).
# Arguments
    - `t`: Time value or vector [s]
    - `rb`: Borehole radius [m]
    - `xy`: Borehole coordinates (nb × 2) [m], e.g. from [`borefield`](@ref)
    - `m`: Ground model (e.g. `FLSModel(150, 4, 3.0, 2e6)`)
    - `nseg`: Number of segments per borehole
# Output
    - `g`: Response array [°C·m/W]
# Reference
    - Dusseault, B., Pasquier, P., & Marcotte, D. (2018). A block matrix formulation for efficient
        g-function construction. Renewable Energy, 121, 249–260.
        https://doi.org/10.1016/j.renene.2017.12.092
    - Cimmino, M., & Bernier, M. (2014). A semi-analytical method to generate g-functions for
        geothermal bore fields. International Journal of Heat and Mass Transfer, 70, 641–650.
        https://doi.org/10.1016/j.ijheatmasstransfer.2013.11.037
    - Cimmino, M. (2015). The effects of borehole thermal resistances and fluid flow rate on the
        g-functions of geothermal bore fields. International Journal of Heat and Mass Transfer,
        91, 1119–1127. https://doi.org/10.1016/j.ijheatmasstransfer.2015.08.041
"""
function gmatrix(t, rb, xy, m::ILSModel)
    r = size(xy, 1) == 1 ? rb : borefield_geometry(xy, rb)[1]
    return ils(t, r, m.ks, m.Cs)
end
function gmatrix(t, rb, xy, m::ICSModel)
    r = size(xy, 1) == 1 ? rb : borefield_geometry(xy, rb)[1]
    return ics(t, r, rb, m.ks, m.Cs)
end
function gmatrix(t, rb, xy, m::FLSModel)
    r = size(xy, 1) == 1 ? rb : borefield_geometry(xy, rb)[1]
    return fls(t, r, m.H, m.D, m.ks, m.Cs)
end
function gmatrix(t, rb, xy, m::MILSModel)
    rb == m.rb || throw(ArgumentError(
        "borehole radius rb ($rb) must match the model's rb ($(m.rb))"))
    size(xy, 1) == 1 && return mils(t, m.rb, 0.0, m.rb, m.ks, m.Cs, m.Cf, m.vD)
    r, θ = borefield_geometry(xy, m.rb)
    return mils(t, r, θ, m.rb, m.ks, m.Cs, m.Cf, m.vD)
end
function gmatrix(t, rb, xy, m::MFLSModel)
    rb == m.rb || throw(ArgumentError(
        "borehole radius rb ($rb) must match the model's rb ($(m.rb))"))
    size(xy, 1) == 1 && return mfls(t, m.rb, 0.0, m.H, m.rb, m.D, m.ks, m.Cs, m.Cf, m.vD)
    r, θ = borefield_geometry(xy, m.rb)
    return mfls(t, r, θ, m.H, m.rb, m.D, m.ks, m.Cs, m.Cf, m.vD)
end
function gmatrix(t, rb, xy::AbstractMatrix{<:Real}, m::FLSModel, nseg::Integer)
    nseg ≥ 1 || throw(ArgumentError("nseg must be ≥ 1"))
    nb = size(xy, 1)
    H, D, ks, Cs = m.H, m.D, m.ks, m.Cs

    Hs   = H / nseg                                   # Segment length (equal-length segments)
    Dseg = [D + (k - 1) * Hs for k in 1:nseg]          # Buried depth of each segment top [m]
    nS   = nb * nseg

    # Horizontal borehole-to-borehole distances
    r, = borefield_geometry(xy, rb)
    rmat = max.(r, float(rb))

    # Precompute the segment-to-segment g-functions for every unique borehole distance
    ru   = unique(round.(vec(rmat); digits = 10))
    ridx = Dict(ru[k] => k for k in eachindex(ru))
    nt   = length(t)
    cache = Array{Vector{Float64}}(undef, length(ru), nseg, nseg)
    for k in eachindex(ru)
        for is in 1:nseg, js in 1:is
            g = fls(t, ru[k], Hs, Dseg[js], Hs, Dseg[is], ks, Cs)
            cache[k, js, is] = g
            cache[k, is, js] = g
        end
    end

    # Segment index → (borehole, segment)
    seg = [(bh, k) for bh in 1:nb for k in 1:nseg]

    g3D = zeros(nt, nS, nS)
    for q in 1:nS
        jb, js = seg[q]                                # emitting segment
        for p in 1:nS
            ib, is = seg[p]                            # receiving segment
            k = ridx[round(rmat[ib, jb]; digits = 10)]
            @views g3D[:, p, q] .= cache[k, js, is]
        end
    end
    return g3D
end
