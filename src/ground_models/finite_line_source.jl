using SpecialFunctions: erf
using QuadGK: quadgk

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
    # Calculate terms
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
    # The lower limit of the integral depends on time
    α = ks / Cs
    lower_lim = inv(sqrt(4 * α * t))
    
    # Perform numerical integration
    integral, _ = quadgk(s -> _fls_integrand(s, H, r, D), lower_lim, T(Inf), rtol = T(1e-6))
    
    return integral / (4 * T(π) * ks)
end

"""
    fls(t, H, r, D, ks, Cs)

Computes the finite line source (FLS) model based on Claesson and Javed (2011). The output is a 
g-function that requires a heat load per unit of borehole length [W/m] to provide the borehole
wall temperature.
# Arguments
    - `t`: Time value or vector [s]
    - `H`: Borehole depth [m]
    - `r`: Radius at which to computed (typically the borehole radius) [m]
        - If `r` is a vector, the output will be a matrix of g-function with columns corresponding
            to different radius and rows to different time steps.
        - If `r` is a 2D array, the output will be a 3D array of g-function with dimensions
            corresponding to (time, x and y) coordinates of the borefield. The radius matrix can
            be computed with `borefield_radius()` from utils.jl.
    - `D`: Buried depth [m]
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
# Output
    - `g`: A g-function corresponding to the borehole wall temperature of the borehole [°Cm/W]
# Reference
    - Claesson, J., & Javed, S. (2011). An analytical method to calculate borehole fluid 
        temperatures for time-scales from minutes to decades. ASHRAE Transactions, 117(PART 2), 
        279–288.
# Example
    g = fls(60:60:3600, 150, 0.076, 4.0, 3.0, 2e6)
"""
function fls(t::Real, H::Real, r::Real, D::Real, ks::Real, Cs::Real)
    # Method for 1 time step and 1 radius
    T = float(promote_type(typeof(t), typeof(H), typeof(r), typeof(D), typeof(ks), typeof(Cs)))
    return _fls(T(t), T(H), T(r), T(D), T(ks), T(Cs))
end
function fls(t::AbstractVector{<:Real}, H::Real, r::Real, D::Real, ks::Real, Cs::Real)
    # Method for multiple time steps and 1 radius.
    # Check type
    T = float(promote_type(eltype(t), typeof(H), typeof(r), typeof(D), typeof(ks), typeof(Cs)))
    t_T  = convert(Vector{T}, t)
    
    # Preallocate and FLS
    g = similar(t_T)
    @inbounds @simd for i in eachindex(t_T)
        g[i] = _fls(t_T[i], T(H), T(r), T(D), T(ks), T(Cs))
    end
    return g
end
function fls(t::Real, H::Real, r::AbstractVector{<:Real}, D::Real, ks::Real, Cs::Real)
    # Method for 1 time step and multiple radius.
    T = float(promote_type(typeof(t), typeof(H), eltype(r), typeof(D), typeof(ks), typeof(Cs)))
    g = similar(r)
    @inbounds @simd for i in eachindex(r)
        g[i] = _fls(T(t), T(H), T(r[i]), T(D), T(ks), T(Cs))
    end
    return g
end
function fls(t::AbstractVector{<:Real}, H::Real, r::AbstractVector{<:Real}, D::Real, ks::Real,
    Cs::Real)
    # Method for multiple time steps and multiple radius in a 2D array. Columns correspond to
    # different radius and rows to different time steps.
    T = float(promote_type(eltype(t), typeof(H), eltype(r), typeof(D), typeof(ks), typeof(Cs)))
    t_T  = convert(Vector{T}, t)
    g = Matrix{T}(undef, length(t_T), length(r))
    @inbounds @simd for i in eachindex(r)
        for j in eachindex(t_T)
            g[j, i] = _fls(t_T[j], T(H), T(r[i]), T(D), T(ks), T(Cs))
        end
    end
    return g
end
function fls(t::AbstractVector{<:Real}, H::Real, r::AbstractArray{<:Real}, D::Real, ks::Real,
    Cs::Real)
    # Method for multiple time steps and multiple radius in a 3D array.
    # The dimensions of the 3D array are (`t`, `x`, `y`).
    # Check type
    T = float(promote_type(eltype(t), typeof(H), eltype(r), typeof(D), typeof(ks), typeof(Cs)))
    # Convert time vector to the right type
    t_T  = convert(Vector{T}, t)
    nt = length(t_T)                # Number of element in the time vector
    # Find unique radius
    nb = size(r, 1)                 # Number of boreholes
    rᵥ = reshape(r, nb * nb)        # Vector of the borefield radius (nb x 1) [m]
    rᵤ = unique(rᵥ)                 # Unique values of the borefield radius (nbᵤ x 1) [m]
    rᵢ = indexin(rᵥ, rᵤ)            # Indices of the unique radius values (nb*nb x 1) [m]
    # Compute a 2D matrix of g-function for unique radius
    g2D = fls(t_T, T(H), rᵤ, T(D), T(ks), T(Cs))
    # Fill a 3D matrix of g-function for all radius
    g3D = zeros(nt, nb, nb)
    for i in 1:nt
        g3D[i, :, :] = reshape(g2D[i, rᵢ], (1, nb, nb)) # Fill a 3D matrix of g-functions
    end
    return g3D
end

function fls_old(t::Union{Real, AbstractVector{<:Real}}, ks::Real, Cs::Real, r::Real, H::Real,
    D::Real)
    # Set initial parameters
    const_π = 1 / sqrt(π)
    nt = length(t)              # Number of element in the time vector
    g = zeros(nt)               # Preallocation of the borehole wall temperature
    α = ks / Cs
    lim_int = 1 ./ sqrt.(4 * α * t)

    """
        integrand_fls(s, r, H, D)

    Integrand of the FLS model. Assumes constant heat flux boundary condition.
    """
    function integrand_fls(s::Real, r::Real, H::Real, D::Real)
        """
            ierf(x)

        Inverse "erf" function used in the FLS model.
        """
        function ierf(x::Real)
            return x * erf(x) - const_π * (1 - exp(-x^2))
        end

        return exp(-r^2 * s^2) * (2 * ierf(H * s) + 2 * ierf(H * s + 2 * D * s) -
                ierf(2 * H * s + 2 * D * s) - ierf(2 * D * s)) / (H * s^2)
    end
    # Compute, in a loop, each value of the fls
    for i in 1:nt
        integral, _ = quadgk(s -> integrand_fls(s, r, H, D), lim_int[i], Inf, rtol = 1e-6)
        g[i] = integral
    end
    return g / (4 * π * ks)
end