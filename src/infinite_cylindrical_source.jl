using SpecialFunctions: besselj0, besselj1, bessely0, bessely1
using QuadGK: quadgk

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
    return _ics(T(t), T(r), T(rc), T(ks), T(Cs))
end
function ics(t::AbstractVector{<:Real}, r::Real, rc::Real, ks::Real, Cs::Real)
    # Method for multiple time steps and 1 radius.
    # Check type
    T = float(promote_type(eltype(t), typeof(r), typeof(rc), typeof(ks), typeof(Cs)))
    t_T  = convert(Vector{T}, t)
    
    # Preallocate and ICS
    g = similar(t_T)
    @inbounds @simd for i in eachindex(t_T)
        g[i] = _ics(t_T[i], T(r), T(rc), T(ks), T(Cs))
    end
    return g
end
function ics(t::Real, r::AbstractVector{<:Real}, rc::Real, ks::Real, Cs::Real)
    # Method for 1 time step and multiple radius.
    T = float(promote_type(typeof(t), eltype(r), typeof(rc), typeof(ks), typeof(Cs)))
    g = similar(r)
    @inbounds @simd for i in eachindex(r)
        g[i] = _ics(T(t), T(r[i]), T(rc), T(ks), T(Cs))
    end
    return g
end
function ics(t::AbstractVector{<:Real}, r::AbstractVector{<:Real}, rc::Real, ks::Real, Cs::Real)
    # Method for multiple time steps and multiple radius in a 2D array. Columns correspond to
    # different radius and rows to different time steps.
    T = float(promote_type(eltype(t), eltype(r), typeof(rc), typeof(ks), typeof(Cs)))
    t_T  = convert(Vector{T}, t)
    g = Matrix{T}(undef, length(t_T), length(r))
    @inbounds @simd for i in eachindex(r)
        for j in eachindex(t_T)
            g[j, i] = _ics(t_T[j], T(r[i]), T(rc), T(ks), T(Cs))
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
    # The dimensions of the 3D array are (`t`, `x`, `y`).
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

"""
    _ics_integrand(s, r̃, t̃)

Computes the integrand of the infinite cylindrical source model.
"""
function _ics_integrand(s::T, r̃::T, t̃::T) where {T<:AbstractFloat}
    if s < 1e-12
        return zero(T)
    end
    # Pre-calculate Bessel terms to keep the expression readable
    j0_rs = besselj0(r̃ * s)
    y1_s  = bessely1(s)
    y0_rs = bessely0(r̃ * s)
    j1_s  = besselj1(s)
    
    return ((exp(-s^2 * t̃) - 1) * (j0_rs * y1_s - y0_rs * j1_s)) / (s^2 * (j1_s^2 + y1_s^2))
end

"""
    _ics(t, r, rc, ks, Cs)

Kernel function for the infinite cylindrical source model based on Carslaw and Jaeger (1959). The
response function is based on an impulse of 1 W/m.
"""
function _ics(t::T, r::T, rc::T, ks::T, Cs::T) where {T<:AbstractFloat}
    r̃ = r / rc
    t̃ = (t * ks) / (Cs * rc^2)

    # For r >> rc the cylinder is indistinguishable from a line source: ICS → ILS
    r̃ > 50 && return _ils(t, r, ks, Cs)

    # Upper integration limit
    s_upper = sqrt(T(50) / t̃)

    if r̃ > T(1.5)
        # The integrand oscillates at spatial frequency r̃/π (from besselj0/bessely0 at r̃·s).
        half_period = T(π) / r̃
        n_breaks = min(floor(Int, s_upper / half_period), 500)
        pts = Vector{T}(undef, n_breaks + 2)
        pts[1] = T(1e-8)
        for i in 1:n_breaks
            pts[i + 1] = i * half_period
        end
        pts[end] = s_upper
        integral, _ = quadgk(s -> _ics_integrand(s, r̃, t̃), pts..., rtol = T(1e-6))
    else
        # r ≈ rc: integrand is smooth (Wronskian identity makes oscillations cancel exactly
        # at r̃ = 1). Standard integration with Inf is fast and well-conditioned.
        integral, _ = quadgk(s -> _ics_integrand(s, r̃, t̃), T(1e-8), T(Inf), rtol = T(1e-6))
    end

    return integral / (T(π)^2 * ks)
end
