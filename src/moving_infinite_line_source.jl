using SpecialFunctions: besseli, besselk, expint
using LinearAlgebra

"""
    mils(t, xy, rb, ks, Cs, Cf, vD)

Compute the moving infinite line source (MILS) model based on Pasquier and Lamarche (2022). The
output is a g-function that requires a heat load per unit of borehole length [W/m] to provide the
borehole wall temperature. The groundwater flow is in the positive x-direction.
# Arguments
    - `t`: Time value or vector [s]
    - `xy`: Coordinates at which to compute [m] — direction-dependent form (Eq. 1).
        - Pass a 2-element vector `[x, y]` for a single point.
        - Pass a `nb×2` matrix of borehole coordinates for a borefield; output is then a
            `(nb × nb)` or `(nt × nb × nb)` g-matrix.
    - `rb`: Borehole radius [m] — required with `xy` to identify inside/outside the borehole.
        - If `r = sqrt(x² + y²) ≤ rb`, uses borehole-wall series scaled by I₀ (Eq. 13 analog).
        - If `r > rb`, uses directional form: ḡ(r) × exp(x·vT/(2α)) / I₀(r·vT/(2α)) (Eq. 1).
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
    - `Cf`: Groundwater volumetric specific heat [J/m³K]
    - `vD`: Darcy velocity (groundwater speed) [m/s]
        - Must not be zero; set to a low value for impervious ground (1e-12) or use `ils`.
# Output
    - `g`: A g-function corresponding to the borehole wall temperature [°Cm/W]
# Reference
    - Pasquier, P., & Lamarche, L. (2022). Analytic expressions for the moving infinite line source
        model. Geothermics, 103, 102413. https://doi.org/10.1016/j.geothermics.2022.102413
# Example
    mils(60:60:3600, [5.0, 0.0], 0.076, 3.0, 2e6, 4.2e6, 1e-6)  # Directional at (5 m, 0 m)
"""
function mils(t::Real, xy::AbstractVector{<:Real}, rb::Real, ks::Real, Cs::Real, Cf::Real,
    vD::Real)
    T = float(promote_type(typeof(t), eltype(xy), typeof(rb), typeof(ks), typeof(Cs), typeof(Cf),
        typeof(vD)))
    return _mils(T(t), T(ks), T(Cs), T(Cf), T(xy[1]), T(xy[2]), T(rb), T(vD))
end
function mils(t::AbstractVector{<:Real}, xy::AbstractVector{<:Real}, rb::Real, ks::Real,
    Cs::Real, Cf::Real, vD::Real)
    T = float(promote_type(eltype(t), eltype(xy), typeof(rb), typeof(ks), typeof(Cs), typeof(Cf),
        typeof(vD)))
    t_T = convert(Vector{T}, t)
    g = similar(t_T)
    @inbounds @simd for i in eachindex(t_T)
        g[i] = _mils(t_T[i], T(ks), T(Cs), T(Cf), T(xy[1]), T(xy[2]), T(rb), T(vD))
    end
    return g
end
function mils(t::Real, xy::AbstractMatrix{<:Real}, rb::Real, ks::Real, Cs::Real, Cf::Real,
    vD::Real)
    T = float(promote_type(typeof(t), eltype(xy), typeof(rb), typeof(ks), typeof(Cs), typeof(Cf),
        typeof(vD)))
    nb = size(xy, 1)
    g2D = Matrix{T}(undef, nb, nb)
    for j in 1:nb
        for i in 1:nb
            dx = T(xy[i, 1] - xy[j, 1])
            dy = T(xy[i, 2] - xy[j, 2])
            g2D[i, j] = _mils(T(t), T(ks), T(Cs), T(Cf), dx, dy, T(rb), T(vD))
        end
    end
    return g2D
end
function mils(t::AbstractVector{<:Real}, xy::AbstractMatrix{<:Real}, rb::Real, ks::Real,
    Cs::Real, Cf::Real, vD::Real)
    T = float(promote_type(eltype(t), eltype(xy), typeof(rb), typeof(ks), typeof(Cs), typeof(Cf),
        typeof(vD)))
    t_T = convert(Vector{T}, t)
    nb = size(xy, 1)
    nt = length(t_T)
    g3D = zeros(T, nt, nb, nb)
    for j in 1:nb
        for i in 1:nb
            dx = T(xy[i, 1] - xy[j, 1])
            dy = T(xy[i, 2] - xy[j, 2])
            for k in eachindex(t_T)
                g3D[k, i, j] = _mils(t_T[k], T(ks), T(Cs), T(Cf), dx, dy, T(rb), T(vD))
            end
        end
    end
    return g3D
end

"""
    _mils_series(t, ks, Cs, Cf, r, vD)

Kernel for the azimuthal mean MILS response (Eq. 3, Pasquier & Lamarche 2022). Returns the
circumferential-average g-function at radius `r` from the line source. Used directly by the
scalar-`r` public overloads, and as the series engine inside `_mils`.
"""
function _mils_series(t::T, ks::T, Cs::T, Cf::T, r::T, vD::T) where {T<:AbstractFloat}
    ns = 5                             # Number of summand used. Below 20 to avoid BigInt issues.
    b = (r * vD * Cf / (4 * ks))^2      # b = (Peclet/4)^2
    τ = (4 * ks / Cs) * (t / r^2)       # τ = 4 * Fo, and Fo = αt/r²

    x = 2 * sqrt(b)
    I0 = besseli(zero(T), x)
    a0 = I0 / (4 * T(π) * ks)

    if τ <= (1 / b)
        # Eq. 20 of Pasquier et Lamarche (2022) (Early Time)
        S1, S2 = zero(T), zero(T)
        for i in 0:(ns-1)
            n = i + 1
            S2 = b^n / factorial(n)^2
            S1 += (-τ)^(i + 1) * factorial(i) * S2
        end
        return a0 * (expint(inv(τ)) * I0 + exp(-inv(τ)) * S1)
    else
        # Eq. 25 of Pasquier et Lamarche (2022) (Late Time)
        S1 = zero(T)
        for i in 0:(ns-1)
            S2 = zero(T)
            for n in 1:ns
                fact = i + n > 18 ? factorial(big(i + n)) : factorial(i + n)
                S2 += b^(n - 1) / T(fact)^2
            end
            S1 += float((factorial(i) / (-τ)^(i + 1)) * S2)
        end
        return a0 * (2 * besselk(zero(T), x) - expint(b * τ) * I0 - exp(-b * τ) * S1)
    end
end

"""
    _mils(t, ks, Cs, Cf, x, y, rb, vD)

Unified direction-dependent kernel for the MILS model (Pasquier & Lamarche 2022). Mirrors the
structure of `_mfls`: switches formula based on whether the evaluation point is inside or outside
the borehole cylinder.
- `r = sqrt(x² + y²) ≤ rb`: uses borehole-wall series `_mils_series(rb)`, scaled by the ratio
  I₀(r·vT/(2α)) / I₀(rb·vT/(2α)) — analogous to Eq. (13) of Guo et al. (2020).
- `r > rb`: directional form — ḡ(r) × exp(x·vT/(2α)) / I₀(r·vT/(2α)) [Eq. 1].
"""
function _mils(t::T, ks::T, Cs::T, Cf::T, x::T, y::T, rb::T, vD::T) where {T<:AbstractFloat}
    α = ks / Cs
    vT = vD * Cf / Cs
    r = sqrt(x^2 + y^2)

    if r ≤ rb
        # Inside borehole: evaluate series at rb, weight by I₀ ratio
        g_rb  = _mils_series(t, ks, Cs, Cf, rb, vD)
        I0_rb = besseli(zero(T), rb * vT / (T(2) * α))
        I0_r  = besseli(zero(T), r  * vT / (T(2) * α))
        return g_rb * I0_r / I0_rb
    else
        # Outside borehole: direction-dependent form (Eq. 1)
        g_r  = _mils_series(t, ks, Cs, Cf, r, vD)
        I0_r = besseli(zero(T), r * vT / (T(2) * α))
        return g_r * exp(x * vT / (T(2) * α)) / I0_r
    end
end
