using SpecialFunctions: besseli, besselk, expint
using QuadGK: quadgk
using LinearAlgebra

"""
    mils(t, r, θ, rb, ks, Cs, Cf, vD)

Compute the moving infinite line source (MILS) model based on Pasquier and Lamarche (2022). The
output is a g-function that requires a heat load per unit of borehole length [W/m] to provide the
borehole wall temperature. The groundwater flow is in the positive x-direction.

Like the isotropic models (`ils`, `fls`), `mils` dispatches on geometry rather than coordinates:
the direction-dependent response is a function of the separation `r` and the flow-relative angle
`θ` only (the point is `x = r·cosθ`, `y = r·sinθ`, and only `x` and `r` enter the model). Pass
scalars for a single point, or matching `nb×nb` matrices, built with [`borefield_geometry`](@ref),
for a borefield.
# Arguments
    - `t`: Time value or vector [s]
    - `r`: Separation distance [m]. Scalar, or an `nb×nb` matrix for a borefield.
        - If `r ≤ rb`, uses the borehole-wall series scaled by I₀ (Eq. 13 analog); the borefield
          self entries (diagonal `r = rb`) therefore give the borehole-wall self response.
        - If `r > rb`, uses the directional form ḡ(r) × exp(x·vT/(2α)) / I₀(r·vT/(2α)) (Eq. 1).
    - `θ`: Flow-relative angle in **degrees**, `θ ∈ [0, 180]`. Scalar, or an `nb×nb` matrix
        matching `r`.
    - `rb`: Borehole radius [m] — identifies inside/outside the borehole.
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
    - `Cf`: Groundwater volumetric specific heat [J/m³K]
    - `vD`: Darcy velocity (groundwater speed) [m/s]
        - Must not be zero, set to a low value for impervious ground (1e-12).
# Output
    - `g`: A g-function corresponding to the borehole wall temperature [°Cm/W]
# Reference
    - Pasquier, P., & Lamarche, L. (2022). Analytic expressions for the moving infinite line source
        model. Geothermics, 103, 102413. https://doi.org/10.1016/j.geothermics.2022.102413
"""
function mils(t::Real, r::Real, θ::Real, rb::Real, ks::Real, Cs::Real, Cf::Real, vD::Real)
    T = float(promote_type(typeof(t), typeof(r), typeof(θ), typeof(rb), typeof(ks), typeof(Cs),
        typeof(Cf), typeof(vD)))
    return _mils(T(t), T(ks), T(Cs), T(Cf), T(r), T(r) * cosd(T(θ)), T(rb), T(vD))
end
function mils(t::AbstractVector{<:Real}, r::Real, θ::Real, rb::Real, ks::Real, Cs::Real,
    Cf::Real, vD::Real)
    T = float(promote_type(eltype(t), typeof(r), typeof(θ), typeof(rb), typeof(ks), typeof(Cs),
        typeof(Cf), typeof(vD)))
    t_T = convert(Vector{T}, t)
    x = T(r) * cosd(T(θ))
    g = similar(t_T)
    @inbounds @simd for i in eachindex(t_T)
        g[i] = _mils(t_T[i], T(ks), T(Cs), T(Cf), T(r), x, T(rb), T(vD))
    end
    return g
end
function mils(t::Real, r::AbstractMatrix{<:Real}, θ::AbstractMatrix{<:Real}, rb::Real, ks::Real,
    Cs::Real, Cf::Real, vD::Real)
    return dropdims(mils([t], r, θ, rb, ks, Cs, Cf, vD); dims = 1)
end
function mils(t::AbstractVector{<:Real}, r::AbstractMatrix{<:Real}, θ::AbstractMatrix{<:Real},
    rb::Real, ks::Real, Cs::Real, Cf::Real, vD::Real)
    @assert size(r) == size(θ) "r and θ must have the same shape"
    T = float(promote_type(eltype(t), eltype(r), eltype(θ), typeof(rb), typeof(ks), typeof(Cs),
        typeof(Cf), typeof(vD)))
    t_T = convert(Vector{T}, t)
    nt = length(t_T)
    nb = size(r, 1)
    pairs = collect(zip(vec(r), vec(θ)))
    uniq = unique(pairs)
    ui = indexin(pairs, uniq)
    gU = Matrix{T}(undef, nt, length(uniq))
    for (u, (rr, tt)) in enumerate(uniq)
        R = T(rr)
        X = R * cosd(T(tt))
        for k in 1:nt
            gU[k, u] = _mils(t_T[k], T(ks), T(Cs), T(Cf), R, X, T(rb), T(vD))
        end
    end
    g3D = zeros(T, nt, nb, nb)
    for k in 1:nt
        @views g3D[k, :, :] .= reshape(gU[k, ui], nb, nb)
    end
    return g3D
end

"""
    _mils(t, ks, Cs, Cf, r, x, rb, vD)

Kernel function for the MILS model (Pasquier & Lamarche 2022). It evaluates the
azimuthal-mean Hantush well function series (Eqs. 20 & 25) and applies the direction-dependent
weighting in one pass, switching formula based on whether the evaluation point is inside or
outside the borehole cylinder (`x = r·cosθ` is the downstream coordinate):
- `r ≤ rb`: the series is evaluated at the borehole wall `rb` and scaled by the ratio
  I₀(r·vT/(2α)) / I₀(rb·vT/(2α)) — analogous to Eq. (13) of Guo et al. (2020). At the self entry
  `r = rb` this ratio is 1, giving the circumferential-mean borehole-wall response.
- `r > rb`: directional (θ-dependent) form — ḡ(r) × exp(x·vT/(2α)) / I₀(r·vT/(2α)) [Eq. 1], so
  downstream points (x > 0) are warmer than upstream ones (x < 0).
"""
function _mils(t::T, ks::T, Cs::T, Cf::T, r::T, x::T, rb::T, vD::T) where {T<:AbstractFloat}
    α  = ks / Cs
    vT = vD * Cf / Cs

    # Radius `rs` at which the azimuthal-mean series is evaluated
    if r ≤ rb
        rs = rb
        weight = besseli(zero(T), r  * vT / (T(2) * α)) /
                 besseli(zero(T), rb * vT / (T(2) * α))
    else
        rs = r
        weight = exp(x * vT / (T(2) * α)) /
                 besseli(zero(T), rs * vT / (T(2) * α))
    end

    # Azimuthal-mean (circumferential-average) response at radius `rs`
    ns = 5                              # Number of summands. Below 20 to avoid BigInt issues.
    b  = (rs * vD * Cf / (4 * ks))^2    # b = (Péclet/4)²
    τ  = (4 * ks / Cs) * (t / rs^2)     # τ = 4·Fo, with Fo = αt/rs²
    xb = 2 * sqrt(b)                    # argument of the Bessel functions
    I0 = besseli(zero(T), xb)
    a0 = I0 / (4 * T(π) * ks)

    if b > one(T)
        W, _ = quadgk(ψ -> exp(-ψ - b / ψ) / ψ, inv(τ), T(Inf), rtol = T(1e-8))
        mean = a0 * W
    elseif τ <= (1 / b)
        # Eq. 20 of Pasquier & Lamarche (2022) (early time).
        S1 = zero(T)
        for m in 0:(ns-1)
            inner = zero(T)
            for n in (m+1):ns
                inner += b^n / factorial(n)^2
            end
            S1 += (-τ)^(m + 1) * factorial(m) * inner
        end
        mean = a0 * (expint(inv(τ)) * I0 + exp(-inv(τ)) * S1)
    else
        # Eq. 25 of Pasquier & Lamarche (2022) (late time): for b ≤ 1 reached only when τ ≥ 1.
        S1 = zero(T)
        for i in 0:(ns-1)
            S2 = zero(T)
            for n in 1:ns
                fact = i + n > 18 ? factorial(big(i + n)) : factorial(i + n)
                S2  += b^(n - 1) / T(fact)^2
            end
            S1 += float((factorial(i) / (-τ)^(i + 1)) * S2)
        end
        mean = a0 * (2 * besselk(zero(T), xb) - expint(b * τ) * I0 - exp(-b * τ) * S1)
    end

    return mean * weight
end
