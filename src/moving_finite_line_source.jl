using SpecialFunctions: erf, besseli
using QuadGK: quadgk

"""
    mfls(t, xy, H, rb, D, ks, Cs, Cf, vD)

Compute the moving finite line source (MFLS) model of Guo et al. (2020), which integrates the
buried depth and groundwater flow (with direction). The output is a g-function that requires a
heat load per unit of borehole length [W/m] to provide the borehole wall temperature. The
groundwater flow is considered to be on the positive x-axis.
# Arguments
    - `t`: Time value or vector [s]
    - `xy`: Coordinates at which to compute [m]
        - Pass a 2-element vector `[x, y]` for a single point.
        - Pass a `nb×2` matrix of borehole coordinates for a borefield; output is then a
            `(nb × nb)` or `(nt × nb × nb)` g-matrix.
        - If `r = sqrt(x² + y²) ≤ rb`, uses borehole-wall integrand with I₀ factor (Eq. 13).
        - If `r > rb`, uses directional form with exp factor (Eq. 10).
    - `H`: Borehole depth [m]
    - `rb`: Borehole radius [m]
    - `D`: Buried depth [m]
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
"""
function mfls(t::Real, xy::AbstractVector{<:Real}, H::Real, rb::Real, D::Real, ks::Real,
    Cs::Real, Cf::Real, vD::Real)
    T = float(promote_type(typeof(t), eltype(xy), typeof(H), typeof(rb), typeof(D), typeof(ks),
        typeof(Cs), typeof(Cf), typeof(vD)))
    return _mfls(T(t), T(H), T(rb), T(D), T(xy[1]), T(xy[2]), T(ks), T(Cs), T(Cf), T(vD))
end
function mfls(t::AbstractVector{<:Real}, xy::AbstractVector{<:Real}, H::Real, rb::Real,
    D::Real, ks::Real, Cs::Real, Cf::Real, vD::Real)
    T = float(promote_type(eltype(t), eltype(xy), typeof(H), typeof(rb), typeof(D), typeof(ks),
        typeof(Cs), typeof(Cf), typeof(vD)))
    t_T = convert(Vector{T}, t)
    g = similar(t_T)
    @inbounds @simd for i in eachindex(t_T)
        g[i] = _mfls(t_T[i], T(H), T(rb), T(D), T(xy[1]), T(xy[2]), T(ks), T(Cs), T(Cf), T(vD))
    end
    return g
end
function mfls(t::Real, xy::AbstractMatrix{<:Real}, H::Real, rb::Real, D::Real, ks::Real,
    Cs::Real, Cf::Real, vD::Real)
    T = float(promote_type(typeof(t), eltype(xy), typeof(H), typeof(rb), typeof(D), typeof(ks),
        typeof(Cs), typeof(Cf), typeof(vD)))
    nb = size(xy, 1)
    g2D = Matrix{T}(undef, nb, nb)
    for j in 1:nb
        for i in 1:nb
            dx = T(xy[i, 1] - xy[j, 1])
            dy = T(xy[i, 2] - xy[j, 2])
            g2D[i, j] = _mfls(T(t), T(H), T(rb), T(D), dx, dy, T(ks), T(Cs), T(Cf), T(vD))
        end
    end
    return g2D
end
function mfls(t::AbstractVector{<:Real}, xy::AbstractMatrix{<:Real}, H::Real, rb::Real,
    D::Real, ks::Real, Cs::Real, Cf::Real, vD::Real)
    T = float(promote_type(eltype(t), eltype(xy), typeof(H), typeof(rb), typeof(D), typeof(ks),
        typeof(Cs), typeof(Cf), typeof(vD)))
    t_T = convert(Vector{T}, t)
    nb = size(xy, 1)
    nt = length(t_T)
    g3D = zeros(T, nt, nb, nb)
    for j in 1:nb
        for i in 1:nb
            dx = T(xy[i, 1] - xy[j, 1])
            dy = T(xy[i, 2] - xy[j, 2])
            for k in eachindex(t_T)
                g3D[k, i, j] = _mfls(t_T[k], T(H), T(rb), T(D), dx, dy, T(ks), T(Cs), T(Cf), T(vD))
            end
        end
    end
    return g3D
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

Unified kernel for the moving finite line source model (Guo et al. 2020). Switches between the
direction-dependent form (Eq. 10) and the circumferential-average form (Eq. 13) based on whether
`r = sqrt(x²+y²)` is outside or inside the borehole radius `rb`.
"""
function _mfls(t::T, H::T, rb::T, D::T, x::T, y::T, ks::T, Cs::T, Cf::T, vD::T
    ) where {T<:AbstractFloat}
    α = ks / Cs
    U = vD * Cf / Cs
    r = sqrt(x^2 + y^2)
    lower_lim = inv(sqrt(T(4) * α * t))

    if r ≤ rb
        # Eq. 13 of Guo et al. (2020) — circumferential average at borehole wall
        r̃ = rb
        I = besseli(zero(T), r * U * inv(T(2) * α))
    else
        # Eq. 10 of Guo et al. (2020) — directional response outside borehole
        r̃ = r
        θ = atan(y, x)
        I = exp(r * U * cos(θ) * inv(T(2) * α))
    end

    integral, _ = quadgk(s -> _mfls_integrand(s, H, r̃, D, α, U), lower_lim, T(Inf), rtol = T(1e-6))

    return (integral * I) / (T(4) * T(π) * ks)
end
