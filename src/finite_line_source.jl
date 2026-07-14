using SpecialFunctions: erf
using QuadGK: quadgk

"""
    fls(t, r, H, D, ks, Cs)

Computes the finite line source (FLS) model based on Claesson and Javed (2011). The output is a
g-function that requires a heat load per unit of borehole length [W/m] to provide the borehole
wall temperature.
# Arguments
    - `t`: Time value or vector [s]
    - `r`: Radius at which to compute [m]
        - If `r` is a vector, the output will be a matrix of g-function with columns corresponding
            to different radii and rows to different time steps.
        - If `r` is a 2D array, the output will be a 3D array of g-function with dimensions
            corresponding to (time, x, y) coordinates of the borefield. The radius matrix can
            be computed with `borefield_geometry()` from borefield.jl.
    - `H`: Borehole depth [m]
    - `D`: Buried depth [m]
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
# Output
    - `g`: A g-function corresponding to the borehole wall temperature of the borehole [°Cm/W]
# Reference
    - Claesson, J., & Javed, S. (2011). An analytical method to calculate borehole fluid
        temperatures for time-scales from minutes to decades. ASHRAE Transactions, 117(PART 2),
        279–288.
"""
function fls(t::Real, r::Real, H::Real, D::Real, ks::Real, Cs::Real)
    T = float(promote_type(typeof(t), typeof(r), typeof(H), typeof(D), typeof(ks), typeof(Cs)))
    return _fls(T(t), T(H), T(r), T(D), T(ks), T(Cs))
end
function fls(t::AbstractVector{<:Real}, r::Real, H::Real, D::Real, ks::Real, Cs::Real)
    T = float(promote_type(eltype(t), typeof(r), typeof(H), typeof(D), typeof(ks), typeof(Cs)))
    t_T = convert(Vector{T}, t)
    g = similar(t_T)
    @inbounds @simd for i in eachindex(t_T)
        g[i] = _fls(t_T[i], T(H), T(r), T(D), T(ks), T(Cs))
    end
    return g
end
function fls(t::Real, r::AbstractVector{<:Real}, H::Real, D::Real, ks::Real, Cs::Real)
    T = float(promote_type(typeof(t), eltype(r), typeof(H), typeof(D), typeof(ks), typeof(Cs)))
    g = similar(r)
    @inbounds @simd for i in eachindex(r)
        g[i] = _fls(T(t), T(H), T(r[i]), T(D), T(ks), T(Cs))
    end
    return g
end
function fls(t::AbstractVector{<:Real}, r::AbstractVector{<:Real}, H::Real, D::Real, ks::Real,
    Cs::Real)
    T = float(promote_type(eltype(t), eltype(r), typeof(H), typeof(D), typeof(ks), typeof(Cs)))
    t_T = convert(Vector{T}, t)
    g = Matrix{T}(undef, length(t_T), length(r))
    @inbounds @simd for i in eachindex(r)
        for j in eachindex(t_T)
            g[j, i] = _fls(t_T[j], T(H), T(r[i]), T(D), T(ks), T(Cs))
        end
    end
    return g
end
function fls(t::Real, r::AbstractArray{<:Real}, H::Real, D::Real, ks::Real,
    Cs::Real)
    T = float(promote_type(typeof(t), eltype(r), typeof(H), typeof(D), typeof(ks), typeof(Cs)))
    nb = size(r, 1)
    rᵥ = reshape(r, nb * nb)
    rᵤ = unique(rᵥ)
    rᵢ = indexin(rᵥ, rᵤ)
    g = fls(T(t), rᵤ, T(H), T(D), T(ks), T(Cs))
    return reshape(g[rᵢ], nb, nb)
end

function fls(t::AbstractVector{<:Real}, r::AbstractArray{<:Real}, H::Real, D::Real, ks::Real,
    Cs::Real)
    # Method for multiple time steps and a 2D radius matrix — builds a 3D g-array (nt × nb × nb).
    T = float(promote_type(eltype(t), eltype(r), typeof(H), typeof(D), typeof(ks), typeof(Cs)))
    t_T = convert(Vector{T}, t)
    nt = length(t_T)                # Number of element in the time vector
    nb = size(r, 1)                 # Number of boreholes
    rᵥ = reshape(r, nb * nb)        # Vector of the borefield radius (nb x 1) [m]
    rᵤ = unique(rᵥ)                 # Unique values of the borefield radius (nbᵤ x 1) [m]
    rᵢ = indexin(rᵥ, rᵤ)            # Indices of the unique radius values (nb*nb x 1) [m]
    g2D = fls(t_T, rᵤ, T(H), T(D), T(ks), T(Cs))
    g3D = zeros(nt, nb, nb)
    for i in 1:nt
        g3D[i, :, :] = reshape(g2D[i, rᵢ], (1, nb, nb))
    end
    return g3D
end

"""
    _ierf(x)

Inverse "erf" function used in the finite line source model.
"""
function _ierf(x::T) where {T<:AbstractFloat}
    return x * erf(x) - inv(sqrt(T(π))) * (one(T) - exp(-x^2))
end

"""
    _fls_integrand(s, H, r, D)

Computes the integrand of the finite line source model.
"""
function _fls_integrand(s::T, H::T, r::T, D::T) where {T<:AbstractFloat}
    fun = (2 * _ierf(H * s)) + (2 * _ierf((H * s) + (2 * D * s))) -
        _ierf((2 * H * s) + (2 * D * s)) - _ierf(2 * D * s)
    return (exp(-r^2 * s^2) * fun) / (H * s^2)
end

"""
    _fls(t, H, r, D, ks, Cs)

Kernel function for the finite line source model based on Claesson and Javed (2011). The response
function is based on an impulse of 1 W/m.
"""
function _fls(t::T, H::T, r::T, D::T, ks::T, Cs::T) where {T<:AbstractFloat}
    α = ks / Cs
    lower_lim = inv(sqrt(4 * α * t))
    integral, _ = quadgk(s -> _fls_integrand(s, H, r, D), lower_lim, T(Inf), rtol = T(1e-6))
    return integral / (4 * T(π) * ks)
end

"""
    fls(t, r, H1, D1, H2, D2, ks, Cs)

Segment-to-segment finite line source (Cimmino & Bernier, 2014): the mean g-function over a
receiving segment `(H2, D2)` produced by a 1 W/m impulse on an emitting segment `(H1, D1)`, at
horizontal separation `r`. This is the general form used to build the segment response matrix of the
BC-III (segment) spatial superposition; the single-argument `fls(t, r, H, D, ks, Cs)` is the special
case `H1 = H2 = H`, `D1 = D2 = D`.
# Arguments
    - `t`: Time value or vector [s]
    - `r`: Horizontal distance between the two segments [m]
    - `H1`, `D1`: Length and buried depth of the emitting segment [m]
    - `H2`, `D2`: Length and buried depth of the receiving segment [m]
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
# Output
    - `g`: g-function averaged over the receiving segment [°Cm/W]
# Reference
    - Cimmino, M., & Bernier, M. (2014). A semi-analytical method to generate g-functions for
        geothermal bore fields. International Journal of Heat and Mass Transfer, 70, 641–650.
        https://doi.org/10.1016/j.ijheatmasstransfer.2013.11.037
"""
function fls(t::Real, r::Real, H1::Real, D1::Real, H2::Real, D2::Real, ks::Real, Cs::Real)
    T = float(promote_type(typeof(t), typeof(r), typeof(H1), typeof(D1), typeof(H2), typeof(D2),
        typeof(ks), typeof(Cs)))
    return _fls_seg(T(t), T(H1), T(D1), T(H2), T(D2), T(r), T(ks), T(Cs))
end

function fls(t::AbstractVector{<:Real}, r::Real, H1::Real, D1::Real, H2::Real, D2::Real, ks::Real,
    Cs::Real)
    T = float(promote_type(eltype(t), typeof(r), typeof(H1), typeof(D1), typeof(H2), typeof(D2),
        typeof(ks), typeof(Cs)))
    t_T = convert(Vector{T}, t)
    g = similar(t_T)
    @inbounds @simd for i in eachindex(t_T)
        g[i] = _fls_seg(t_T[i], T(H1), T(D1), T(H2), T(D2), T(r), T(ks), T(Cs))
    end
    return g
end

"""
    _fls_seg_integrand(s, H1, D1, H2, D2, r)

Integrand of the segment-to-segment finite line source solution. `(H1, D1)` are the length and
buried depth of the emitting segment, `(H2, D2)` those of the receiving segment, and `r` the
horizontal distance between the two segments.
"""
function _fls_seg_integrand(s::T, H1::T, D1::T, H2::T, D2::T, r::T) where {T<:AbstractFloat}
    Ireal = _ierf((D2 - D1 + H2) * s) - _ierf((D2 - D1) * s) +
        _ierf((D2 - D1 - H1) * s) - _ierf((D2 - D1 + H2 - H1) * s)
    Iimag = _ierf((D2 + D1 + H2) * s) - _ierf((D2 + D1) * s) +
        _ierf((D2 + D1 + H1) * s) - _ierf((D2 + D1 + H2 + H1) * s)
    return (exp(-r^2 * s^2) * (Ireal + Iimag)) / (H2 * s^2)
end

"""
    _fls_seg(t, H1, D1, H2, D2, r, ks, Cs)

Kernel for the segment-to-segment finite line source. Returns the mean g-function over the receiving
segment `(H2, D2)` due to a 1 W/m impulse on the emitting segment `(H1, D1)`, at horizontal
separation `r`.
"""
function _fls_seg(t::T, H1::T, D1::T, H2::T, D2::T, r::T, ks::T, Cs::T) where {T<:AbstractFloat}
    α = ks / Cs
    lower_lim = inv(sqrt(4 * α * t))
    integral, _ = quadgk(s -> _fls_seg_integrand(s, H1, D1, H2, D2, r), lower_lim, T(Inf),
        rtol = T(1e-6), atol = T(1e-8))
    return integral / (4 * T(π) * ks)
end
