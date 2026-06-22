using SpecialFunctions: expinti

"""
    ils(t, r, ks, Cs)

Compute the infinite line source (ILS) model based on Ingersol (1954). The output is a g-function
that requires a heat load per unit of borehole length [W/m] to provide the borehole wall
temperature.
# Arguments
    - `t`: Time value or vector [s]
    - `r`: Radius at which to computed (typically the borehole radius) [m]
        - If `r` is a vector, the output will be a matrix of g-function with columns corresponding
            to different radius and rows to different time steps.
        - If `r` is a 2D array, the output will be a 3D array of g-function with dimensions
            corresponding to (time, x and y) coordinates of the borefield. The radius matrix can
            be computed with `borefield_radius()` from borefield.jl.
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
# Output
    - `g`: A g-function corresponding to the borehole wall temperature of the borehole [°Cm/W]
# Reference
    - Ingersol, L. R. (1948). Theory of the ground pipe heat source for the heat pump. 
        ASHVE Journal Section, Heating, Piping and Air Conditioning.
# Example
    g = ils(60:60:3600, 3.0, 2e6, 0.076)
"""
function ils(t::Real, r::Real, ks::Real, Cs::Real)
    # Method for 1 time step and 1 radius
    T = float(promote_type(typeof(t), typeof(r), typeof(ks), typeof(Cs)))
    return _ils(T(t), T(r), T(ks), T(Cs))
end
function ils(t::AbstractVector{<:Real}, r::Real, ks::Real, Cs::Real)
    # Method for multiple time steps and 1 radius.
    # Check type
    T = float(promote_type(eltype(t), typeof(r), typeof(ks), typeof(Cs)))
    t_T  = convert(Vector{T}, t)

    # Preallocate and ILS
    g = similar(t_T)
    @inbounds @simd for i in eachindex(t_T)
        g[i] = _ils(t_T[i], T(r), T(ks), T(Cs))
    end
    return g
end
function ils(t::Real, r::AbstractVector{<:Real}, ks::Real, Cs::Real)
    # Method for 1 time step and multiple radius.
    T = float(promote_type(typeof(t), eltype(r), typeof(ks), typeof(Cs)))
    g = similar(r)
    @inbounds @simd for i in eachindex(r)
        g[i] = _ils(T(t), T(r[i]), T(ks), T(Cs))
    end
    return g
end
function ils(t::AbstractVector{<:Real}, r::AbstractVector{<:Real}, ks::Real, Cs::Real)
    # Method for multiple time steps and multiple radius in a 2D array. Columns correspond to
    # different radius and rows to different time steps.
    T = float(promote_type(eltype(t), eltype(r), typeof(ks), typeof(Cs)))
    t_T  = convert(Vector{T}, t)
    g = Matrix{T}(undef, length(t_T), length(r))
    @inbounds @simd for i in eachindex(r)
        for j in eachindex(t_T)
            g[j, i] = _ils(t_T[j], T(r[i]), T(ks), T(Cs))
        end
    end
    return g
end
function ils(t::Real, r::AbstractArray{<:Real}, ks::Real, Cs::Real)
    T = float(promote_type(typeof(t), eltype(r), typeof(ks), typeof(Cs)))
    nb = size(r, 1)
    rᵥ = reshape(r, nb * nb)
    rᵤ = unique(rᵥ)
    rᵢ = indexin(rᵥ, rᵤ)
    g1D = ils(T(t), rᵤ, T(ks), T(Cs))
    return reshape(g1D[rᵢ], nb, nb)
end
function ils(t::AbstractVector{<:Real}, r::AbstractArray{<:Real}, ks::Real, Cs::Real)
    # Method for multiple time steps and multiple radius in a 3D array.
    T = float(promote_type(eltype(t), eltype(r), typeof(ks), typeof(Cs)))
    t_T  = convert(Vector{T}, t)
    nt = length(t_T)                # Number of element in the time vector
    nb = size(r, 1)                 # Number of boreholes
    rᵥ = reshape(r, nb * nb)        # Vector of the borefield radius (nb x 1) [m]
    rᵤ = unique(rᵥ)                 # Unique values of the borefield radius (nbᵤ x 1) [m]
    rᵢ = indexin(rᵥ, rᵤ)            # Indices of the unique radius values (nb*nb x 1) [m]
    g2D = ils(t_T, rᵤ, T(ks), T(Cs))
    g3D = zeros(nt, nb, nb)
    for i in 1:nt
        g3D[i, :, :] = reshape(g2D[i, rᵢ], (1, nb, nb)) # Fill a 3D matrix of g-functions
    end
    return g3D
end

"""
    _ils(t, r, ks, Cs)

Kernel function for the infinite line source model based on Ingersol (1954). The response function 
is based on an impulse of 1 W/m.
"""
@inline function _ils(t::T, r::T, ks::T, Cs::T) where {T<:AbstractFloat}
    α = ks / Cs
    x = -r^2 / (4 * α * t)
    return -expinti(x) / (4 * T(π) * ks)
end