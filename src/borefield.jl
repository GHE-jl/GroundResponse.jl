using LinearAlgebra

"""
    borefield_geometry(xy, rb)

Pairwise geometry of a borefield: the **distance matrix** and the **flow-relative angle matrix**.
These are the two quantities every ground model consumes:
    - isotropic models (ILS, ICS, FLS) use only the distance matrix `r`;
    - moving models (MILS, MFLS) additionally use the angle matrix `θ`, because under groundwater
        flow the pairwise response is asymmetric (`g[i,j] ≠ g[j,i]`) and depends on both the
        separation `rᵢⱼ` and the angle `θᵢⱼ` of the source→receiver direction relative to the
        flow (+x).
Convention: row `i` is the **receiver** and column `j` the **source**, matching `g[k,i,j]`. The
angle is measured for the source→receiver vector, so a receiver directly downstream of the source
has `θ = 0°` (warmest) and one directly upstream has `θ = 180°` (coolest). Angles are returned in
**degrees**.
# Arguments
    - `xy`: Borehole coordinates (nb × 2) [m]. Groundwater flow is assumed along +x; to orient the
        field relative to the flow, generate it with the `angle` keyword of the borefield
        generators (e.g. [`borefield_rectangle`](@ref)).
    - `rb`: Borehole radius [m]. Placed on the diagonal of `r`.
# Outputs
    - `r`: Pairwise distance matrix (nb × nb) [m] (diagonal entries equal `rb`).
    - `θ`: Pairwise angle matrix (nb × nb) [°]; `θ[i,j] = acosd((xᵢ − xⱼ)/r[i,j])` ∈ [0, 180];
        `0` on the diagonal.
# Inspecting the borefield further
To recover the distinct `(distance, angle)` combinations, an index from each borehole pair to its
combination, and the multiplicity of each, work directly from `r` and `θ`:
```julia
r, θ = borefield_geometry(xy, rb)
combos = round.(hcat(vec(r), vec(θ)); digits = 10)        # (r, θ) for every (i, j) pair
keys   = unique(eachrow(combos))                          # distinct (distance, angle) combinations
idx    = reshape(indexin(eachrow(combos), keys), size(r)) # nb×nb: which combination each pair is
counts = [count(==(k), eachrow(combos)) for k in keys]    # how many pairs share each (Σ = nb²)
```
`keys` is what an optimised solver would evaluate once, `idx` scatters those results back onto the
`nb×nb` grid, and a large `counts` entry flags a highly symmetric field (many equivalent pairs).
"""
function borefield_geometry(xy::AbstractMatrix{<:Real}, rb::Real)
    nb = size(xy, 1)
    dx = xy[:, 1] .- xy[:, 1]'                          # dx[i,j] = xᵢ − xⱼ (receiver − source)
    r  = hypot.(dx, xy[:, 2] .- xy[:, 2]')
    θ  = [i == j ? 0.0 : acosd(clamp(dx[i, j] / r[i, j], -1.0, 1.0)) for i in 1:nb, j in 1:nb]
    r  = r + Diagonal(fill(float(rb), nb))              # self-distance on the diagonal is rb
    return r, θ
end

"""
    borefield(shape, args...; angle = 0)

Unified entry point for generating borehole coordinates. `shape` selects the layout; the remaining
arguments and the `angle` keyword are forwarded unchanged to the corresponding layout function.
# Arguments
    | `shape`            | Layout function              | Signature(s)                         |
    |:-------------------|:-----------------------------|:-------------------------------------|
    | `:rectangle`       | `borefield_rectangle`        | `(nx, ny, B)` or `(nx, ny, Bx, By)`  |
    | `:line`            | `borefield_line`             | `(n, B)`                             |
    | `:circle`          | `borefield_circle`           | `(nb, R)`                            |
    | `:L`               | `borefield_L`                | `(n1, n2, B)` or `(n1, n2, B1, B2)`  |
    | `:U`               | `borefield_U`                | `(nx, ny, B)` or `(nx, ny, Bx, By)`  |
    | `:open_rectangle`  | `borefield_open_rectangle`   | `(nx, ny, B)` or `(nx, ny, Bx, By)`  |
# Keywords
    - `angle`: Rotation of the borefield, counter-clockwise about the origin [°] (default `0`)
# Output
    - `xy`: Borehole coordinates (nb × 2) [m]
"""
function borefield(shape::Symbol, args...; angle::Real = 0)
    shape === :rectangle      && return borefield_rectangle(args...; angle)
    shape === :line           && return borefield_line(args...; angle)
    shape === :circle         && return borefield_circle(args...; angle)
    shape === :L              && return borefield_L(args...; angle)
    shape === :U              && return borefield_U(args...; angle)
    shape === :open_rectangle && return borefield_open_rectangle(args...; angle)
    throw(ArgumentError(
        "Unknown borefield shape :$shape. " *
        "Valid shapes: :rectangle, :line, :circle, :L, :U, :open_rectangle"
    ))
end

"""
    borefield_rectangle(nx, ny, B; angle = 0)
    borefield_rectangle(nx, ny, Bx, By; angle = 0)

Generate a rectangular (or square) grid of borehole coordinates.
# Arguments
    - `nx`: Number of boreholes in the x direction [-]
    - `ny`: Number of boreholes in the y direction [-]
    - `B`: Uniform spacing between boreholes [m]
    - `Bx`: Spacing in the x direction [m]
    - `By`: Spacing in the y direction [m]
# Keywords
    - `angle`: Rotation of the borefield, counter-clockwise about the origin [°] (default `0`)
# Output
    - `xy`: Borehole coordinates (nb × 2) [m], `nb = nx * ny`
"""
function borefield_rectangle(nx::Integer, ny::Integer, B::Real; angle::Real = 0)
    xy = hcat([[i, j] for i in 1:nx for j in 1:ny]...)' .* B .- B
    return _rotate(xy, angle)
end
function borefield_rectangle(nx::Integer, ny::Integer, Bx::Real, By::Real; angle::Real = 0)
    xy = hcat([[i, j] for i in 1:nx for j in 1:ny]...)' .* [Bx By] .- [Bx By]
    return _rotate(xy, angle)
end

"""
    borefield_line(n, B; angle = 0)

Generate a single row of `n` boreholes equally spaced along the x-axis.
# Arguments
    - `n`: Number of boreholes [-]
    - `B`: Spacing between boreholes [m]
# Keywords
    - `angle`: Rotation of the borefield, counter-clockwise about the origin [°] (default `0`)
# Output
    - `xy`: Borehole coordinates (n × 2) [m]
"""
function borefield_line(n::Integer, B::Real; angle::Real = 0)
    xy = hcat((0:n-1) .* B, zeros(n))
    return _rotate(xy, angle)
end

"""
    borefield_circle(nb, R; angle = 0)

Generate `nb` boreholes evenly distributed on a circle of radius `R`, starting at angle 0 (positive
x-axis).
# Arguments
    - `nb`: Number of boreholes on the circle [-]
    - `R`: Radius of the circle [m]
# Keywords
    - `angle`: Rotation of the borefield, counter-clockwise about the origin [°] (default `0`)
# Output
    - `xy`: Borehole coordinates (nb × 2) [m]
"""
function borefield_circle(nb::Integer, R::Real; angle::Real = 0)
    θ  = (2π / nb) .* (0:nb-1)
    xy = hcat(R .* cos.(θ), R .* sin.(θ))
    return _rotate(xy, angle)
end

"""
    borefield_L(n1, n2, B; angle = 0)
    borefield_L(n1, n2, B1, B2; angle = 0)

Generate an L-shaped borefield. The horizontal arm has `n1` boreholes along the x-axis and the
vertical arm has `n2` boreholes along the y-axis; they share the corner at the origin. Total number
of boreholes is `n1 + n2 - 1`.
# Arguments
    - `n1`: Number of boreholes in the horizontal arm (including the shared corner) [-]
    - `n2`: Number of boreholes in the vertical arm (including the shared corner) [-]
    - `B`: Uniform spacing [m]
    - `B1`: Spacing along the horizontal arm [m]
    - `B2`: Spacing along the vertical arm [m]
# Keywords
    - `angle`: Rotation of the borefield, counter-clockwise about the origin [°] (default `0`)
# Output
    - `xy`: Borehole coordinates ((n1+n2-1) × 2) [m]
"""
function borefield_L(n1::Integer, n2::Integer, B::Real; angle::Real = 0)
    h_arm = hcat((0:n1-1) .* B, zeros(n1))
    v_arm = hcat(zeros(n2-1), (1:n2-1) .* B)
    return _rotate(vcat(h_arm, v_arm), angle)
end
function borefield_L(n1::Integer, n2::Integer, B1::Real, B2::Real; angle::Real = 0)
    h_arm = hcat((0:n1-1) .* B1, zeros(n1))
    v_arm = hcat(zeros(n2-1), (1:n2-1) .* B2)
    return _rotate(vcat(h_arm, v_arm), angle)
end

"""
    borefield_U(nx, ny, B; angle = 0)
    borefield_U(nx, ny, Bx, By; angle = 0)

Generate a U-shaped borefield: `nx` boreholes across the base and `ny` boreholes up each side
(open at the top). Total number of boreholes is `nx + 2*(ny-1)`.
# Arguments
    - `nx`: Number of boreholes along the base (≥ 2) [-]
    - `ny`: Number of boreholes up each side (including the base corners) [-]
    - `B`: Uniform spacing [m]
    - `Bx`: Spacing along the base [m]
    - `By`: Spacing along the sides [m]
# Keywords
    - `angle`: Rotation of the borefield, counter-clockwise about the origin [°] (default `0`)
# Output
    - `xy`: Borehole coordinates ((nx + 2*(ny-1)) × 2) [m]
"""
function borefield_U(nx::Integer, ny::Integer, B::Real; angle::Real = 0)
    bottom = hcat((0:nx-1) .* B,                  zeros(nx))
    left   = hcat(zeros(ny-1),                     (1:ny-1) .* B)
    right  = hcat(fill(Float64((nx-1)*B), ny-1),   (1:ny-1) .* B)
    return _rotate(vcat(bottom, left, right), angle)
end
function borefield_U(nx::Integer, ny::Integer, Bx::Real, By::Real; angle::Real = 0)
    bottom = hcat((0:nx-1) .* Bx,                  zeros(nx))
    left   = hcat(zeros(ny-1),                      (1:ny-1) .* By)
    right  = hcat(fill(Float64((nx-1)*Bx), ny-1),   (1:ny-1) .* By)
    return _rotate(vcat(bottom, left, right), angle)
end

"""
    borefield_open_rectangle(nx, ny, B; angle = 0)
    borefield_open_rectangle(nx, ny, Bx, By; angle = 0)

Generate a hollow rectangular borefield — boreholes only on the perimeter of an `nx × ny` grid.
Total number of boreholes is `2*(nx + ny - 2)`. Requires `nx ≥ 2` and `ny ≥ 2`.
# Arguments
    - `nx`: Number of boreholes along each horizontal side [-]
    - `ny`: Number of boreholes along each vertical side [-]
    - `B`: Uniform spacing [m]
    - `Bx`: Spacing along horizontal sides [m]
    - `By`: Spacing along vertical sides [m]
# Keywords
    - `angle`: Rotation of the borefield, counter-clockwise about the origin [°] (default `0`)
# Output
    - `xy`: Borehole coordinates (2*(nx+ny-2) × 2) [m]
"""
function borefield_open_rectangle(nx::Integer, ny::Integer, B::Real; angle::Real = 0)
    bottom = hcat((0:nx-1) .* B,                   zeros(nx))
    top    = hcat((0:nx-1) .* B,                   fill(Float64((ny-1)*B), nx))
    left   = hcat(zeros(ny-2),                      (1:ny-2) .* B)
    right  = hcat(fill(Float64((nx-1)*B), ny-2),    (1:ny-2) .* B)
    return _rotate(vcat(bottom, top, left, right), angle)
end
function borefield_open_rectangle(nx::Integer, ny::Integer, Bx::Real,
    By::Real; angle::Real = 0)
    bottom = hcat((0:nx-1) .* Bx,                   zeros(nx))
    top    = hcat((0:nx-1) .* Bx,                   fill(Float64((ny-1)*By), nx))
    left   = hcat(zeros(ny-2),                       (1:ny-2) .* By)
    right  = hcat(fill(Float64((nx-1)*Bx), ny-2),    (1:ny-2) .* By)
    return _rotate(vcat(bottom, top, left, right), angle)
end

"""
    _rotate(xy, angle)

Rotate the borehole coordinates `xy` (nb × 2) [m] counter-clockwise about the origin by `angle`
[°]. Used by the borefield generators, so that a field can be oriented relative to the groundwater flow
(+x) of the moving models, or to the site for a plot.
"""
function _rotate(xy::AbstractMatrix{<:Real}, angle::Real)
    angle == 0 && return xy
    c, s = cosd(angle), sind(angle)
    return hcat(c .* xy[:, 1] .- s .* xy[:, 2], s .* xy[:, 1] .+ c .* xy[:, 2])
end
