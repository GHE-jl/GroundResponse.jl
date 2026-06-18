using SpecialFunctions: besseli, besselk, expint
using LinearAlgebra

"""
    _mils(t, ks, Cs, Cf, r, vD)

Kernel function of the moving infinite line source based on Pasquier et Lamarche (2022). The
response is based on an impulse of 1 W/m.
"""
function _mils(t::T, ks::T, Cs::T, Cf::T, r::T, vD::T) where {T<:AbstractFloat}
    # Initial inputs
    ns = 5                             # Number of summand used. Below 20 to avoid BigInt issues.
    b = (r * vD * Cf / (4 * ks))^2      # b = (Peclet/4)^2
    τ = (4 * ks / Cs) * (t / r^2)       # τ = 4 * Fo, and Fo = αt/r²
    
    # Pre-compute common Bessel and coefficient terms
    x = 2 * sqrt(b)
    I0 = besseli(zero(T), x)
    a0 = I0 / (4 * T(π) * ks)
    
    # MILS
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
                i + n > 18 ? factorial(big(i + n)) : factorial(i + n)
                S2 += b^(n - 1) / factorial(i + n)^2
            end
            S1 += float((factorial(i) / (-τ)^(i + 1)) * S2)
        end
        return a0 * (2 * besselk(zero(T), x) - expint(b * τ) * I0 - exp(-b * τ) * S1)
    end
end

"""
    mils(t, ks, Cs, Cf, r, vD)

Compute the moving infinite line source (MILS) model based on Pasquier et Lamarche (2022). The 
output is a g-function that requires a heat load per unit of borehole length [W/m] to provide the 
borehole wall temperature.
# Arguments
    - `t`: Time value or vector [s]
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
    - `Cf`: Groundwater volumetric specific heat [J/m³K]
    - `r`: Radius at which to computed (typically the borehole radius) [m]
    - `vD`: Darcy velocity (groundwater speed) [m/s]
        - Must not be zero, set to low value for impervious (1e-12) or use the "ils" model.
# Output
    - `g`: A g-function corresponding to the borehole wall temperature of the borehole [°Cm/W]
# Reference
    - Pasquier, P., & Lamarche, L. (2022). Analytic expressions for the moving infinite line source
        model. Geothermics, 103, 102413. https://doi.org/10.1016/j.geothermics.2022.102413
# Example
    mils(60:60:3600, 3.0, 2e6, 4.2e6, 0.076, 1e-6)
"""
function mils(t::Real, ks::Real, Cs::Real, Cf::Real, r::Real, vD::Real)
    T = float(promote_type(typeof(t), typeof(ks), typeof(Cs), typeof(Cf), typeof(r), typeof(vD)))
    return _mils(T(t), T(ks), T(Cs), T(Cf), T(r), T(vD))
end

function mils(t::AbstractVector{<:Real}, ks::Real, Cs::Real, Cf::Real, r::Real, vD::Real)
    # Check type
    T = float(promote_type(eltype(t), typeof(ks), typeof(Cs), typeof(Cf), typeof(r), typeof(vD)))
    t_T  = convert(Vector{T}, t)

    # Preallocate and MILS
    g = similar(t_T)
    @inbounds @simd for i in eachindex(t_T)
        g[i] = _mils(t_T[i], T(ks), T(Cs), T(Cf), T(r), T(vD))
    end
    return g
end

function mils!(g::AbstractVector{T}, t::AbstractVector, ks::Real, Cs::Real, Cf::Real, r::Real,
    vD::Real) where {T<:AbstractFloat}
    # Check for same vector length
    @assert length(g) == length(t)

    # Convert parameters once to the target type
    ks_T, Cs_T, Cf_T, r_T, vD_T = T(ks), T(Cs), T(Cf), T(r), T(vD)

    # MILS
    @inbounds @simd for i in eachindex(t)
        g[i] = _mils(T(t[i]), ks_T, Cs_T, Cf_T, r_T, vD_T)
    end
    return g
end