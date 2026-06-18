using SpecialFunctions: erf, besseli
using QuadGK: quadgk

"""
    _ierf(x)

Inverse "erf" function used in the moving finite line source model.
Equation 12 of Guo et al. (2020).
"""
function _ierf(x::T) where {T<:AbstractFloat}
    return x * erf(x) - inv(sqrt(T(π))) * (one(T) - exp(-x^2))
end

"""
    _mfls_integrand(s, H, r, D, α, U)

Computes the integrand of the finite line source model. Equations 10 and 13 of Guo et al. (2020).
"""
function _mfls_integrand(s::T, H::T, r::T, D::T, α::T, U::T) where {T<:AbstractFloat}
    # Equation 11 of Guo et al. (2020)
    fun = (2 * _ierf(H * s)) + (2 * _ierf((H * s) + (2 * D * s))) - 
        _ierf((2 * H * s) + (2 * D * s)) - _ierf(2 * D * s)
    
    return exp(-U^2 * inv(T(16) * α^2 * s^2) - (r^2 * s^2)) * fun / (H * s^2)
end

"""
    _mfls(t, H, rb, D, x, y, ks, Cs, Cf, vD)

Kernel function for the moving finite line source model based on Guo et al. (2020). The response 
function is based on an impulse of 1 W/m.
"""
function _mfls(t::T, H::T, rb::T, D::T, x::T, y::T, ks::T, Cs::T, Cf::T, vD::T
    ) where {T<:AbstractFloat}
    # Initial parameters
    α = ks / Cs
    U = vD * Cf / Cs
    r = sqrt(x^2 + y^2)
    lower_lim = inv(sqrt(T(4) * α * t))

    # Determine how to compute the integrand based on radius
    if r ≤ rb
        # Eq. 13 of Guo et al. (2020) (Inside the borehole)
        r̃ = rb
        I = besseli(zero(T), r * U * inv(T(2) * α))
    else
        # Eq. 10 of Guo et al. (2020) (Outside the borehole)
        r̃ = r
        θ = atan(y, x)
        I = exp(r * U * cos(θ) * inv(T(2) * α))
    end
    
    # Numerical integration
    integral, _ = quadgk(s -> _mfls_integrand(s, H, r̃, D, α, U), lower_lim, T(Inf), rtol = T(1e-6))
    
    return (integral * I) / (T(4) * T(π) * ks)
end

"""
    mfls(t, H, rb, D, xy, ks, Cs, Cf, vD)

Compute the moving finite line source (MFLS) model of Guo et al. (2020), which integrates the 
buried depth and groundwater flow (with direction). The output is a g-function that requires a
heat load per unit of borehole length [W/m] to provide the borehole wall temperature. The
groundwater flow is considered to be on the positive x-axis.
# Arguments
    - `t`: Time value or vector [s]
    - `H`: Borehole depth [m]
    - `rb`: Borehole radius [m]
    - `D`: Buried depth [m]
    - `xy`: Coordinates at which to compute (2x1) [m]
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
    - `Cf`: Groundwater volumetric specific heat [J/m³K]
    - `vD`: Uniform Darcy velocity [m/s]
        - Must not be zero, set to low value for impervious (1e-12)
# Output
    - `g`: A g-function corresponding to the borehole wall temperature of the borehole [°Cm/W]
# Reference
    - Guo, Y., Hu, X., Banks, J., & Liu, W. V. (2020). Considering buried depth in the moving
        finite line source model for vertical borehole heat exchangers—A new solution. Energy and 
        Buildings, 214, 109859. https://doi.org/10.1016/j.enbuild.2020.109859
# Example
    mfls(60:60:3600, 150, 0.076, 4, [0, 0], 3.0, 2e6, 4.2e6, 1e-6)  # At borehole wall
    mfls(60:60:3600, 150, 0.076, 4, [5, 0], 3.0, 2e6, 4.2e6, 1e-6)  # 5 m from the line source
"""
function mfls(t::Real, H::Real, rb::Real, D::Real, xy::AbstractMatrix{<:Real}, ks::Real, Cs::Real,
    Cf::Real, vD::Real)
    # Method for 1 time step and 1 radius
    T = float(promote_type(typeof(t), typeof(H), typeof(rb), typeof(D), eltype(xy), typeof(ks),
        typeof(Cs), typeof(Cf), typeof(vD)))
    return _mfls(T(t), T(H), T(rb), T(D), T(xy[1]), T(xy[2]), T(ks), T(Cs), T(Cf), T(vD))
end

function mfls(t::AbstractVector{<:Real}, H::Real, rb::Real, D::Real, xy::AbstractMatrix{<:Real},
    ks::Real, Cs::Real, Cf::Real, vD::Real)
    # Method for multiple time steps and 1 xy coordinate.
    # Check type
    T = float(promote_type(eltype(t), typeof(H), typeof(rb), typeof(D), eltype(xy), typeof(ks),
        typeof(Cs), typeof(Cf), typeof(vD)))
    t_T  = convert(Vector{T}, t)
    
    # Preallocate and MFLS
    g = similar(t_T)
    @inbounds @simd for i in eachindex(t_T)
        g[i] = _mfls(t_T[i], T(H), T(rb), T(D), T(xy[1]), T(xy[2]), T(ks), T(Cs), T(Cf), T(vD))
    end
    return g
end
function mfls(t::Real, H::Real, rb::Real, D::Real, xy::AbstractMatrix{<:Real}, ks::Real, Cs::Real,
    Cf::Real, vD::Real)
    # Method for 1 time step and multiple xy coordinates. xy is a 2×n matrix where each column
    # is an [x, y] coordinate pair.
    T = float(promote_type(typeof(t), typeof(H), typeof(rb), typeof(D), eltype(xy), typeof(ks),
        typeof(Cs), typeof(Cf), typeof(vD)))
    n_coords = size(xy, 2)
    g = Vector{T}(undef, n_coords)
    @inbounds for i in 1:n_coords
        g[i] = _mfls(T(t), T(H), T(rb), T(D), T(xy[1, i]), T(xy[2, i]), T(ks), T(Cs), T(Cf), T(vD))
    end
    return g
end
function mfls(t::AbstractVector{<:Real}, H::Real, rb::Real, D::Real, xy::AbstractMatrix{<:Real},
    ks::Real, Cs::Real, Cf::Real, vD::Real)
    # Method for multiple time steps and multiple xy coordinates. xy is a 2×n matrix where each
    # column is an [x, y] coordinate pair.
    T = float(promote_type(eltype(t), typeof(H), typeof(rb), typeof(D), eltype(xy), typeof(ks),
        typeof(Cs), typeof(Cf), typeof(vD)))
    t_T  = convert(Vector{T}, t)
    n_coords = size(xy, 2)
    g = Matrix{T}(undef, length(t_T), n_coords)
    @inbounds for i in 1:n_coords
        for j in eachindex(t_T)
            g[j, i] = _mfls(t_T[j], T(H), T(rb), T(D), T(xy[1, i]), T(xy[2, i]), T(ks), T(Cs), T(Cf), T(vD))
        end
    end
    return g
end