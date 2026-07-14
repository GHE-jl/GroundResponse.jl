using SpecialFunctions: erf, besseli
using QuadGK: quadgk

"""
    mfls(t, r, θ, H, rb, D, ks, Cs, Cf, vD)

Compute the moving finite line source (MFLS) model of Guo et al. (2020), which integrates the
buried depth and groundwater flow (with direction). The output is a g-function that requires a
heat load per unit of borehole length [W/m] to provide the borehole wall temperature. The
groundwater flow is considered to be on the positive x-axis.

Like the isotropic models (`ils`, `fls`), `mfls` dispatches on geometry rather than coordinates:
the response is a function of the separation `r` and the flow-relative angle `θ` only. Pass scalars
for a single point, or matching `nb×nb` matrices — built with [`borefield_geometry`](@ref) — for a
borefield.
# Arguments
    - `t`: Time value or vector [s]
    - `r`: Separation distance [m]. Scalar, or an `nb×nb` matrix for a borefield.
        - If `r ≤ rb`, uses the borehole-wall integrand with I₀ factor (Eq. 13); the borefield
          self entries (diagonal `r = rb`) therefore give the borehole-wall self response.
        - If `r > rb`, uses the directional form with exp factor (Eq. 10).
    - `θ`: Flow-relative angle in **degrees**, `θ ∈ [0, 180]`. Scalar, or an `nb×nb` matrix
        matching `r`. Downstream `θ = 0` is warmest, upstream `θ = 180` coolest.
    - `H`: Borehole depth [m]
    - `rb`: Borehole radius [m]
    - `D`: Buried depth [m]
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
    - `Cf`: Groundwater volumetric specific heat [J/m³K]
    - `vD`: Uniform Darcy velocity [m/s]
        - Must not be zero, set to low value for impervious (1e-12)
# Output
    - `g`: A g-function corresponding to the borehole wall temperature of the borehole [°Cm/W]. A
        borefield input gives an `(nb × nb)` matrix (scalar `t`) or `(nt × nb × nb)` array (vector `t`).
# Reference
    - Guo, Y., Hu, X., Banks, J., & Liu, W. V. (2020). Considering buried depth in the moving
        finite line source model for vertical borehole heat exchangers—A new solution. Energy and
        Buildings, 214, 109859. https://doi.org/10.1016/j.enbuild.2020.109859
"""
function mfls(t::Real, r::Real, θ::Real, H::Real, rb::Real, D::Real, ks::Real,
    Cs::Real, Cf::Real, vD::Real)
    T = float(promote_type(typeof(t), typeof(r), typeof(θ), typeof(H), typeof(rb), typeof(D),
        typeof(ks), typeof(Cs), typeof(Cf), typeof(vD)))
    return _mfls(T(t), T(H), T(rb), T(D), T(r), T(r) * cosd(T(θ)), T(ks), T(Cs), T(Cf), T(vD))
end
function mfls(t::AbstractVector{<:Real}, r::Real, θ::Real, H::Real, rb::Real,
    D::Real, ks::Real, Cs::Real, Cf::Real, vD::Real)
    T = float(promote_type(eltype(t), typeof(r), typeof(θ), typeof(H), typeof(rb), typeof(D),
        typeof(ks), typeof(Cs), typeof(Cf), typeof(vD)))
    t_T = convert(Vector{T}, t)
    x   = T(r) * cosd(T(θ))
    g = similar(t_T)
    @inbounds @simd for i in eachindex(t_T)
        g[i] = _mfls(t_T[i], T(H), T(rb), T(D), T(r), x, T(ks), T(Cs), T(Cf), T(vD))
    end
    return g
end
function mfls(t::Real, r::AbstractMatrix{<:Real}, θ::AbstractMatrix{<:Real}, H::Real, rb::Real,
    D::Real, ks::Real, Cs::Real, Cf::Real, vD::Real)
    return dropdims(mfls([t], r, θ, H, rb, D, ks, Cs, Cf, vD); dims = 1)
end
function mfls(t::AbstractVector{<:Real}, r::AbstractMatrix{<:Real}, θ::AbstractMatrix{<:Real},
    H::Real, rb::Real, D::Real, ks::Real, Cs::Real, Cf::Real, vD::Real)
    @assert size(r) == size(θ) "r and θ must have the same shape"
    T = float(promote_type(eltype(t), eltype(r), eltype(θ), typeof(H), typeof(rb), typeof(D),
        typeof(ks), typeof(Cs), typeof(Cf), typeof(vD)))
    t_T = convert(Vector{T}, t)
    nt  = length(t_T)
    nb  = size(r, 1)

    # Direction-dependent response is asymmetric under advection but shares the same value across
    # borehole pairs with the same (distance, angle). Evaluate each unique (r, θ) response only
    # once, then scatter (Rose et al. 2026 strategy). This is a large saving for the MFLS model,
    # whose kernel performs a quadrature per evaluation. Mirrors fls's unique-radius collapse.
    pairs = collect(zip(vec(r), vec(θ)))
    uniq  = unique(pairs)
    ui    = indexin(pairs, uniq)
    gU    = Matrix{T}(undef, nt, length(uniq))
    for (u, (rr, tt)) in enumerate(uniq)
        R = T(rr)
        X = R * cosd(T(tt))
        for k in 1:nt
            gU[k, u] = _mfls(t_T[k], T(H), T(rb), T(D), R, X, T(ks), T(Cs), T(Cf), T(vD))
        end
    end
    g3D = zeros(T, nt, nb, nb)
    for k in 1:nt
        @views g3D[k, :, :] .= reshape(gU[k, ui], nb, nb)
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
    _mfls(t, H, rb, D, r, x, ks, Cs, Cf, vD)

Unified kernel for the moving finite line source model (Guo et al. 2020). Switches between the
direction-dependent form (Eq. 10) and the circumferential-average form (Eq. 13) based on whether
the separation `r` is outside or inside the borehole radius `rb` (`x = r·cosθ` is the downstream
coordinate). At the self entry `r = rb` the Eq. 13 factor is `I₀(rb·U/2α)`, giving the
circumferential-mean borehole-wall response.
"""
function _mfls(t::T, H::T, rb::T, D::T, r::T, x::T, ks::T, Cs::T, Cf::T, vD::T
    ) where {T<:AbstractFloat}
    α = ks / Cs
    U = vD * Cf / Cs
    lower_lim = inv(sqrt(T(4) * α * t))

    if r ≤ rb
        # Eq. 13 of Guo et al. (2020) — circumferential average at borehole wall
        r̃ = rb
        I = besseli(zero(T), r * U * inv(T(2) * α))
    else
        # Eq. 10 of Guo et al. (2020) — directional response outside borehole (r·cosθ = x)
        r̃ = r
        I = exp(x * U * inv(T(2) * α))
    end

    integral, _ = quadgk(s -> _mfls_integrand(s, H, r̃, D, α, U), lower_lim, T(Inf), rtol = T(1e-6))

    return (integral * I) / (T(4) * T(π) * ks)
end
