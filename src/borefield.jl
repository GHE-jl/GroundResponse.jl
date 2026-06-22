using LinearAlgebra

"""
    borefield_radius(xy, rb)

Compute the pairwise radius matrix, flattened vector, unique values, and indices for a borefield.
Used in spatial superposition to evaluate the ground thermal response at each unique inter-borehole
distance.
# Arguments
    - `xy`: Borehole coordinates (nb × 2) [m]
    - `rb`: Borehole radius [m]
# Outputs
    - `r`: Pairwise radius matrix (nb × nb) [m] — diagonal entries equal `rb`
    - `rᵥ`: Flattened vector of `r` (nb² × 1) [m]
    - `rᵤ`: Unique radius values (nbᵤ × 1) [m]
    - `rᵢ`: Index of each entry of `rᵥ` into `rᵤ` (nb² × 1)
    - `θ`: Azimuth angle of each borehole from the origin (nb,) [rad]
    - `nb`: Number of boreholes [-]
"""
function borefield_radius(xy::AbstractArray{<:Real}, rb::Real)
    r = sqrt.(sum(abs2, xy, dims=2) .+ sum(abs2, xy, dims=2)' .- 2 * (xy * xy'))
    nb = size(xy, 1)
    r  = r + Diagonal(rb * ones(nb))
    rᵥ = reshape(r, nb * nb)
    rᵤ = unique(rᵥ)
    rᵢ = indexin(rᵥ, rᵤ)
    θ  = atan.(xy[:, 2], xy[:, 1])
    return r, rᵥ, rᵤ, rᵢ, θ, nb
end

# ---------------------------------------------------------------------------
# Borefield layout generators
# ---------------------------------------------------------------------------

"""
    borefield_rectangle(nx, ny, B)
    borefield_rectangle(nx, ny, Bx, By)

Generate a rectangular (or square) grid of borehole coordinates.
# Arguments
    - `nx`: Number of boreholes in the x direction [-]
    - `ny`: Number of boreholes in the y direction [-]
    - `B`: Uniform spacing between boreholes [m]
    - `Bx`: Spacing in the x direction [m]
    - `By`: Spacing in the y direction [m]
# Output
    - `xy`: Borehole coordinates (nb × 2) [m], `nb = nx * ny`
"""
function borefield_rectangle(nx::Integer, ny::Integer, B::Real)
    xy = hcat([[i, j] for i in 1:nx for j in 1:ny]...)' .* B .- B
    return xy
end
function borefield_rectangle(nx::Integer, ny::Integer, Bx::Real, By::Real)
    xy = hcat([[i, j] for i in 1:nx for j in 1:ny]...)' .* [Bx By] .- [Bx By]
    return xy
end

"""
    borefield_line(n, B)

Generate a single row of `n` boreholes equally spaced along the x-axis.
# Arguments
    - `n`: Number of boreholes [-]
    - `B`: Spacing between boreholes [m]
# Output
    - `xy`: Borehole coordinates (n × 2) [m]
"""
function borefield_line(n::Integer, B::Real)
    xy = hcat((0:n-1) .* B, zeros(n))
    return xy
end

"""
    borefield_circle(nb, R)

Generate `nb` boreholes evenly distributed on a circle of radius `R`, starting at angle 0 (positive
x-axis).
# Arguments
    - `nb`: Number of boreholes on the circle [-]
    - `R`: Radius of the circle [m]
# Output
    - `xy`: Borehole coordinates (nb × 2) [m]
"""
function borefield_circle(nb::Integer, R::Real)
    θ  = (2π / nb) .* (0:nb-1)
    xy = hcat(R .* cos.(θ), R .* sin.(θ))
    return xy
end

"""
    borefield_L(n1, n2, B)
    borefield_L(n1, n2, B1, B2)

Generate an L-shaped borefield. The horizontal arm has `n1` boreholes along the x-axis and the
vertical arm has `n2` boreholes along the y-axis; they share the corner at the origin. Total number
of boreholes is `n1 + n2 - 1`.
# Arguments
    - `n1`: Number of boreholes in the horizontal arm (including the shared corner) [-]
    - `n2`: Number of boreholes in the vertical arm (including the shared corner) [-]
    - `B`: Uniform spacing [m]
    - `B1`: Spacing along the horizontal arm [m]
    - `B2`: Spacing along the vertical arm [m]
# Output
    - `xy`: Borehole coordinates ((n1+n2-1) × 2) [m]
"""
function borefield_L(n1::Integer, n2::Integer, B::Real)
    h_arm = hcat((0:n1-1) .* B, zeros(n1))
    v_arm = hcat(zeros(n2-1), (1:n2-1) .* B)
    return vcat(h_arm, v_arm)
end
function borefield_L(n1::Integer, n2::Integer, B1::Real, B2::Real)
    h_arm = hcat((0:n1-1) .* B1, zeros(n1))
    v_arm = hcat(zeros(n2-1), (1:n2-1) .* B2)
    return vcat(h_arm, v_arm)
end

"""
    borefield_U(nx, ny, B)
    borefield_U(nx, ny, Bx, By)

Generate a U-shaped borefield: `nx` boreholes across the base and `ny` boreholes up each side
(open at the top). Total number of boreholes is `nx + 2*(ny-1)`.
# Arguments
    - `nx`: Number of boreholes along the base (≥ 2) [-]
    - `ny`: Number of boreholes up each side (including the base corners) [-]
    - `B`: Uniform spacing [m]
    - `Bx`: Spacing along the base [m]
    - `By`: Spacing along the sides [m]
# Output
    - `xy`: Borehole coordinates ((nx + 2*(ny-1)) × 2) [m]
"""
function borefield_U(nx::Integer, ny::Integer, B::Real)
    bottom = hcat((0:nx-1) .* B,                  zeros(nx))
    left   = hcat(zeros(ny-1),                     (1:ny-1) .* B)
    right  = hcat(fill(Float64((nx-1)*B), ny-1),   (1:ny-1) .* B)
    return vcat(bottom, left, right)
end
function borefield_U(nx::Integer, ny::Integer, Bx::Real, By::Real)
    bottom = hcat((0:nx-1) .* Bx,                  zeros(nx))
    left   = hcat(zeros(ny-1),                      (1:ny-1) .* By)
    right  = hcat(fill(Float64((nx-1)*Bx), ny-1),   (1:ny-1) .* By)
    return vcat(bottom, left, right)
end

"""
    borefield_open_rectangle(nx, ny, B)
    borefield_open_rectangle(nx, ny, Bx, By)

Generate a hollow rectangular borefield — boreholes only on the perimeter of an `nx × ny` grid.
Total number of boreholes is `2*(nx + ny - 2)`. Requires `nx ≥ 2` and `ny ≥ 2`.
# Arguments
    - `nx`: Number of boreholes along each horizontal side [-]
    - `ny`: Number of boreholes along each vertical side [-]
    - `B`: Uniform spacing [m]
    - `Bx`: Spacing along horizontal sides [m]
    - `By`: Spacing along vertical sides [m]
# Output
    - `xy`: Borehole coordinates (2*(nx+ny-2) × 2) [m]
"""
function borefield_open_rectangle(nx::Integer, ny::Integer, B::Real)
    bottom = hcat((0:nx-1) .* B,                   zeros(nx))
    top    = hcat((0:nx-1) .* B,                   fill(Float64((ny-1)*B), nx))
    left   = hcat(zeros(ny-2),                      (1:ny-2) .* B)
    right  = hcat(fill(Float64((nx-1)*B), ny-2),    (1:ny-2) .* B)
    return vcat(bottom, top, left, right)
end
function borefield_open_rectangle(nx::Integer, ny::Integer, Bx::Real, By::Real)
    bottom = hcat((0:nx-1) .* Bx,                   zeros(nx))
    top    = hcat((0:nx-1) .* Bx,                   fill(Float64((ny-1)*By), nx))
    left   = hcat(zeros(ny-2),                       (1:ny-2) .* By)
    right  = hcat(fill(Float64((nx-1)*Bx), ny-2),    (1:ny-2) .* By)
    return vcat(bottom, top, left, right)
end

# ---------------------------------------------------------------------------
# Unified entry point
# ---------------------------------------------------------------------------

"""
    borefield(shape, args...)

Unified entry point for generating borehole coordinates. `shape` selects the layout; the remaining
arguments are forwarded unchanged to the corresponding layout function.

| `shape`            | Layout function              | Signature(s)                         |
|:-------------------|:-----------------------------|:-------------------------------------|
| `:rectangle`       | `borefield_rectangle`        | `(nx, ny, B)` or `(nx, ny, Bx, By)` |
| `:line`            | `borefield_line`             | `(n, B)`                             |
| `:circle`          | `borefield_circle`           | `(nb, R)`                            |
| `:L`               | `borefield_L`                | `(n1, n2, B)` or `(n1, n2, B1, B2)` |
| `:U`               | `borefield_U`                | `(nx, ny, B)` or `(nx, ny, Bx, By)` |
| `:open_rectangle`  | `borefield_open_rectangle`   | `(nx, ny, B)` or `(nx, ny, Bx, By)` |

# Example
```julia
xy = borefield(:rectangle, 5, 5, 5.0)   # 5×5 grid, 5 m spacing
xy = borefield(:circle, 8, 10.0)        # 8 boreholes on a circle of radius 10 m
xy = borefield(:L, 4, 3, 5.0)           # L-shape, 4 boreholes × 3 boreholes
```
"""
function borefield(shape::Symbol, args...)
    shape === :rectangle      && return borefield_rectangle(args...)
    shape === :line           && return borefield_line(args...)
    shape === :circle         && return borefield_circle(args...)
    shape === :L              && return borefield_L(args...)
    shape === :U              && return borefield_U(args...)
    shape === :open_rectangle && return borefield_open_rectangle(args...)
    throw(ArgumentError(
        "Unknown borefield shape :$shape. " *
        "Valid shapes: :rectangle, :line, :circle, :L, :U, :open_rectangle"
    ))
end
