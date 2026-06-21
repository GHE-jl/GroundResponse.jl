using LinearAlgebra

"""
    borefield_radius(xy, rb)

Function that computes a radius matrix, vector, unique values and indices of a borefield given
the coordinates of each borehole. This helps to compute the g-function of a borefield using spatial
superposition, as it allows to compute the ground thermal response at each unique radius of the
borefield.
# Arguments
    - `xy`: Matrix of borehole coordinates where the line source is at (0,0) (nb x 2) [m]
        - E.g.: [0 0] (to have a matrix input).
    - `rb`: Borehole radius [m]
# Outputs
    - `r`: Matrix of the borefield radius (nb x nb) [m]
    - `rᵥ`: Vector of the borefield radius ((nb x nb) x 1) [m]
    - `rᵤ`: Unique values of the borefield radius (nbᵤ x 1) [m]
    - `rᵢ`: Indices of the unique radius values (nb*nb x 1) [m]
    - `θ`: Angle of the boreholes from the origin (0,0) (nb,) [rad]
    - `nb`: Number of boreholes [-]
"""
function borefield_radius(xy::AbstractArray{<:Real}, rb::Real)
    # Analyse the radius of the borefield
    r = sqrt.(sum(abs2, xy, dims=2) .+ sum(abs2, xy, dims=2)' .- 2 * (xy * xy'))

    # Outputs
    nb = size(xy, 1)                # Number of boreholes
    r = r + Diagonal(rb * ones(nb)) # Matrix of the borefield radius (nb x nb) [m]
    rᵥ = reshape(r, nb * nb)        # Vector of the borefield radius (nb x 1) [m]
    rᵤ = unique(rᵥ)                 # Unique values of the borefield radius (nbᵤ x 1) [m]
    rᵢ = indexin(rᵥ, rᵤ)            # Indices of the unique radius values (nb*nb x 1) [m]
    θ = atan.(xy[:, 2], xy[:, 1])   # Angle of each borehole from the origin (nb,) [rad]
    return r, rᵥ, rᵤ, rᵢ, θ, nb
end

"""
    borefield_xy(nx, ny, B)
    borefield_xy(nx, ny, Bx, By)

Function that generates the coordinates of a rectangular borefield given the number of boreholes in
the `x` and `y` directions and the spacing between them. For a square borefield, `Bx` and `By` are
equal.
# Arguments
    - `nx`: Number of boreholes in the x direction [-]
    - `ny`: Number of boreholes in the y direction [-]
    - `B`: Spacing between boreholes in both directions [m]
    - `Bx`: Spacing between boreholes in the x direction [m]
    - `By`: Spacing between boreholes in the y direction [m]
# Output
    - `xy`: Matrix of borehole coordinates (nb x 2) [m]
"""
function borefield_xy(nx::Integer, ny::Integer, B::Real)
    xy = hcat([[i, j] for i in 1:nx for j in 1:ny]...)' .* B .- B
    return xy
end
function borefield_xy(nx::Integer, ny::Integer, Bx::Real, By::Real)
    xy = hcat([[i, j] for i in 1:nx for j in 1:ny]...)' .* [Bx By] .- [Bx By]
    return xy
end