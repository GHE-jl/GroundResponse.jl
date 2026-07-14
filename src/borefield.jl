using LinearAlgebra

"""
    borefield_geometry(xy, rb)

Pairwise geometry of a borefield: the distance matrix, the flow-relative angle matrix, and the
**unique (distance, angle) combinations**. This helper feeds every ground model:
  - isotropic models (ILS, ICS, FLS) consume only the distance matrix `r`;
  - moving models (MILS, MFLS) additionally consume the angle matrix `θ`, because under
    groundwater flow the pairwise response is asymmetric (`g[i,j] ≠ g[j,i]`) and depends on both
    the separation `rᵢⱼ` and the angle `θᵢⱼ` of the source→receiver direction relative to the
    flow (+x).
Convention: row `i` is the **receiver** and column `j` the **source**, matching `g[k,i,j]`. The
angle is measured for the source→receiver vector, so a receiver directly downstream of the source
has `θ = 0°` (warmest) and one directly upstream has `θ = 180°` (coolest). Angles are returned in
**degrees**.
# Arguments
    - `xy`: Borehole coordinates (nb × 2) [m]. Groundwater flows along is assumed along +x.
    - `rb`: Borehole radius [m]. Placed on the diagonal of `r`.
# Outputs
    - `r`: Pairwise distance matrix (nb × nb) [m] (diagonal entries equal `rb`).
    - `θ`: Pairwise angle matrix (nb × nb) [°]; `θ[i,j] = acosd((xᵢ − xⱼ)/r[i,j])` ∈ [0, 180];
        `0` on the diagonal.
    - `keys`: Unique `(r, θ)` pairs (nu × 2) [m, °], including the diagonal `(rb, 0)`.
    - `idx`: Index matrix (nb × nb) mapping each `(i,j)` to its row in `keys`. Lets a caller
        evaluate each unique response once and scatter it.
"""
function borefield_geometry(xy::AbstractArray{<:Real}, rb::Real)
    nb = size(xy, 1)
    r  = sqrt.(max.(sum(abs2, xy, dims=2) .+ sum(abs2, xy, dims=2)' .- 2 * (xy * xy'), 0.0))
    dx = xy[:, 1] .- xy[:, 1]'                          # dx[i,j] = xᵢ − xⱼ (receiver − source)
    θ  = [i == j ? 0.0 : acosd(clamp(dx[i, j] / r[i, j], -1.0, 1.0)) for i in 1:nb, j in 1:nb]
    r  = r + Diagonal(fill(float(rb), nb))              # self-distance on the diagonal → rb

    keys   = Tuple{Float64,Float64}[]
    keymap = Dict{Tuple{Float64,Float64},Int}()
    idx    = zeros(Int, nb, nb)
    for j in 1:nb, i in 1:nb
        key = (round(r[i, j]; digits=10), round(θ[i, j]; digits=10))
        idx[i, j] = get!(keymap, key) do
            push!(keys, key)
            length(keys)
        end
    end
    keymat = isempty(keys) ? zeros(0, 2) : reduce(vcat, ([k[1] k[2]] for k in keys))
    return r, θ, keymat, idx
end

"""
    borefield(shape, args...)

Unified entry point for generating borehole coordinates. `shape` selects the layout; the remaining
arguments are forwarded unchanged to the corresponding layout function.
# Arguments
| `shape`            | Layout function              | Signature(s)                         |
|:-------------------|:-----------------------------|:-------------------------------------|
| `:rectangle`       | `borefield_rectangle`        | `(nx, ny, B)` or `(nx, ny, Bx, By)`  |
| `:line`            | `borefield_line`             | `(n, B)`                             |
| `:circle`          | `borefield_circle`           | `(nb, R)`                            |
| `:L`               | `borefield_L`                | `(n1, n2, B)` or `(n1, n2, B1, B2)`  |
| `:U`               | `borefield_U`                | `(nx, ny, B)` or `(nx, ny, Bx, By)`  |
| `:open_rectangle`  | `borefield_open_rectangle`   | `(nx, ny, B)` or `(nx, ny, Bx, By)`  |
# Output
    - `xy`: Borehole coordinates (nb × 2) [m]
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
