using SpecialFunctions: besselj0, besselj1, bessely0, bessely1

"""
    ics(t, r, rc, ks, Cs)

Computes the infinite cylindrical source (ICS) model based on Carslaw and Jaeger (1959). The output
is a g-function that requires a heat load per unit of borehole length [W/m] to provide the borehole
wall temperature.
# Arguments
    - `t`: Time value or vector [s]
    - `r`: Radius at which to compute (typically the borehole radius) [m]
        - If `r` is a vector, the output will be a matrix of g-function with columns corresponding
            to different radius and rows to different time steps.
        - If `r` is a 2D array, the output will be a 3D array of g-function with dimensions
            corresponding to (time, x and y) coordinates of the borefield. The radius matrix can
            be computed with `borefield_geometry()` from borefield.jl.
    - `rc`: Radius of the cylinder [m]
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
# Output
    - `g`: A g-function corresponding to the borehole wall temperature of the borehole [°Cm/W]
# Reference
    - Carslaw, H. S., & Jaeger, J. C. (1959). Conduction of Heat in Solids (2nd ed.).
        Oxford: Clarendon Press.
"""
function ics(t::Real, r::Real, rc::Real, ks::Real, Cs::Real)
    # Method for 1 time step and 1 radius
    T = float(promote_type(typeof(t), typeof(r), typeof(rc), typeof(ks), typeof(Cs)))
    grid = _ics_grid(T)
    return _ics(T(t), T(r), T(rc), T(ks), T(Cs), grid)
end
function ics(t::AbstractVector{<:Real}, r::Real, rc::Real, ks::Real, Cs::Real)
    # Method for multiple time steps and 1 radius.
    # Check type
    T = float(promote_type(eltype(t), typeof(r), typeof(rc), typeof(ks), typeof(Cs)))
    t_T  = convert(Vector{T}, t)
    
    # Preallocate and ICS
    grid = _ics_grid(T)
    g = similar(t_T)
    @inbounds for i in eachindex(t_T)
        g[i] = _ics(t_T[i], T(r), T(rc), T(ks), T(Cs), grid)
    end
    return g
end
function ics(t::Real, r::AbstractVector{<:Real}, rc::Real, ks::Real, Cs::Real)
    # Method for 1 time step and multiple radius.
    T = float(promote_type(typeof(t), eltype(r), typeof(rc), typeof(ks), typeof(Cs)))
    grid = _ics_grid(T)
    g = similar(r, T)
    @inbounds for i in eachindex(r)
        g[i] = _ics(T(t), T(r[i]), T(rc), T(ks), T(Cs), grid)
    end
    return g
end
function ics(t::AbstractVector{<:Real}, r::AbstractVector{<:Real}, rc::Real, ks::Real, Cs::Real)
    # Method for multiple time steps and multiple radius in a 2D array. Columns correspond to
    # different radius and rows to different time steps.
    T = float(promote_type(eltype(t), eltype(r), typeof(rc), typeof(ks), typeof(Cs)))
    t_T  = convert(Vector{T}, t)
    grid = _ics_grid(T)
    g = Matrix{T}(undef, length(t_T), length(r))
    @inbounds for i in eachindex(r)
        for j in eachindex(t_T)
            g[j, i] = _ics(t_T[j], T(r[i]), T(rc), T(ks), T(Cs), grid)
        end
    end
    return g
end
function ics(t::Real, r::AbstractArray{<:Real}, rc::Real, ks::Real, Cs::Real)
    T = float(promote_type(typeof(t), eltype(r), typeof(rc), typeof(ks), typeof(Cs)))
    nb = size(r, 1)
    rᵥ = reshape(r, nb * nb)
    rᵤ = unique(rᵥ)
    rᵢ = indexin(rᵥ, rᵤ)
    g1D = ics(T(t), rᵤ, T(rc), T(ks), T(Cs))
    return reshape(g1D[rᵢ], nb, nb)
end
function ics(t::AbstractVector{<:Real}, r::AbstractArray{<:Real}, rc::Real, ks::Real, Cs::Real)
    # Method for multiple time steps and multiple radius in a 3D array.
    # Check type
    T = float(promote_type(eltype(t), eltype(r), typeof(rc), typeof(ks), typeof(Cs)))
    t_T  = convert(Vector{T}, t)
    nt = length(t_T)                # Number of element in the time vector
    nb = size(r, 1)                 # Number of boreholes
    rᵥ = reshape(r, nb * nb)        # Vector of the borefield radius (nb x 1) [m]
    rᵤ = unique(rᵥ)                 # Unique values of the borefield radius (nbᵤ x 1) [m]
    rᵢ = indexin(rᵥ, rᵤ)            # Indices of the unique radius values (nb*nb x 1) [m]
    g2D = ics(t_T, rᵤ, T(rc), T(ks), T(Cs))
    g3D = zeros(nt, nb, nb)
    for i in 1:nt
        g3D[i, :, :] = reshape(g2D[i, rᵢ], (1, nb, nb)) # Fill a 3D matrix of g-functions
    end
    return g3D
end

# Radius ratio r̃ = r/rc above which the cylinder is indistinguishable from a line source and the
# ICS integral is replaced by the (closed-form, constant-time) ILS. This both removes the most
# expensive integrals and matches the physics: ICS → ILS in the far field.
const _ICS_ILS_RATIO = 20.0

# Number of nodes in the fixed logarithmic quadrature grid.
const _ICS_GRID_N = 100_000

"""
    _ICSGrid{T}

Precomputed, source-independent data for the fixed-grid quadrature of the ICS integral. The
integration nodes `B` and every factor that depends only on `B` (not on the radius, time or
cylinder radius) are computed once and reused for every `(t, r)` evaluation. This gives the model a
constant evaluation cost regardless of the radius ratio `r̃ = r/rc`, unlike an adaptive quadrature
whose cost explodes with the oscillation frequency of the Bessel terms (∝ `r̃`).
- `B`  : integration nodes (log-spaced from 1e-10 to 5e3).
- `BJ` : `besselj1.(B)`.
- `BY` : `bessely1.(B)`.
- `w`  : trapezoidal weight divided by the `B`-only denominator, i.e.
         `dB / (B^2 * (BJ^2 + BY^2))`.
"""
struct _ICSGrid{T<:AbstractFloat}
    B::Vector{T}
    BJ::Vector{T}
    BY::Vector{T}
    w::Vector{T}
end

"""
    _ics_grid(T, n)

Build the fixed logarithmic quadrature grid used by [`_ics`](@ref).
"""
function _ics_grid(::Type{T}, n::Int = _ICS_GRID_N) where {T<:AbstractFloat}
    B = exp10.(range(T(-10), log10(T(5e3)), length = n))
    BJ = besselj1.(B)
    BY = bessely1.(B)

    # Trapezoidal weights on the non-uniform mesh (matches SCI.m).
    w = similar(B)
    @inbounds begin
        w[1]   = (B[2] - B[1]) / 2
        w[end] = (B[end] - B[end-1]) / 2
        for i in 2:n-1
            w[i] = (B[i+1] - B[i-1]) / 2
        end
        # Fold in the B-only denominator so it is not recomputed per evaluation.
        for i in 1:n
            w[i] /= B[i]^2 * (BJ[i]^2 + BY[i]^2)
        end
    end
    return _ICSGrid{T}(B, BJ, BY, w)
end

"""
    _ics(t, r, rc, ks, Cs, grid)

Kernel function for the infinite cylindrical source model based on Carslaw and Jaeger (1959). The
response function is based on an impulse of 1 W/m. The integral is evaluated on the fixed
logarithmic `grid` (see [`_ics_grid`](@ref)), giving a constant cost independent of `r/rc`. For
`r/rc > $(_ICS_ILS_RATIO)` the cylinder is indistinguishable from a line source and the closed-form
ILS is returned instead.
"""
function _ics(t::T, r::T, rc::T, ks::T, Cs::T, grid::_ICSGrid{T}) where {T<:AbstractFloat}
    r̃ = r / rc

    # For r >> rc the cylinder is indistinguishable from a line source: ICS → ILS.
    r̃ > T(_ICS_ILS_RATIO) && return _ils(t, r, ks, Cs)

    t̃ = (t * ks) / (Cs * rc^2)

    B  = grid.B
    BJ = grid.BJ
    BY = grid.BY
    w  = grid.w

    acc = zero(T)
    @inbounds @simd for i in eachindex(B)
        b = B[i]
        num = (exp(-b^2 * t̃) - 1) * (besselj0(r̃ * b) * BY[i] - bessely0(r̃ * b) * BJ[i])
        term = num * w[i]
        acc += ifelse(isfinite(term), term, zero(T))
    end

    return acc / (T(π)^2 * ks)
end
