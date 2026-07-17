# Borefields

Spatial superposition needs the coordinates of every borehole in the field. This page documents the
layout generators that produce those coordinates and the [`borefield_geometry`](@ref) helper that
turns them into the pairwise distance and angle matrices consumed by [`successive_flux`](@ref) and
[`bloc_matrix`](@ref).

## Coordinate convention

Every layout function returns an ``n_b \times 2`` matrix of ``[x \ y]`` coordinates in metres. That
matrix is exactly what [`ground_response`](@ref) and the superposition methods expect:

```julia
xy = borefield(:rectangle, 3, 4, 6.0)    # 3×4 grid, 6 m spacing → 12×2 matrix
g  = ground_response(t, rb, xy, model)   # superposition applied automatically
```

A single borehole is just `[0.0 0.0]` (a `1×2` matrix), which routes `ground_response` to the
direct single-borehole branch.

## The unified entry point

[`borefield`](@ref) forwards its arguments to the matching layout generator according to the first
`Symbol` argument. The individual generators are also exported, so the two calls below are
equivalent:

```julia
xy = borefield(:circle, 8, 10.0)
xy = borefield_circle(8, 10.0)
```

## Available layouts

| Call | Layout | Boreholes |
|---|---|---|
| `borefield(:rectangle, nx, ny, B)` | filled `nx×ny` grid, uniform spacing `B` | `nx·ny` |
| `borefield(:rectangle, nx, ny, Bx, By)` | filled grid, independent x/y spacing | `nx·ny` |
| `borefield(:line, n, B)` | single row along the x-axis | `n` |
| `borefield(:circle, nb, R)` | evenly spaced on a circle of radius `R` | `nb` |
| `borefield(:L, n1, n2, B)` | L-shape, shared corner at the origin | `n1+n2−1` |
| `borefield(:L, n1, n2, B1, B2)` | L-shape, independent arm spacings | `n1+n2−1` |
| `borefield(:U, nx, ny, B)` | U-shape: base of `nx` + `ny` up each side | `nx+2(ny−1)` |
| `borefield(:U, nx, ny, Bx, By)` | U-shape, independent spacings | `nx+2(ny−1)` |
| `borefield(:open_rectangle, nx, ny, B)` | hollow rectangle, perimeter only | `2(nx+ny−2)` |
| `borefield(:open_rectangle, nx, ny, Bx, By)` | hollow rectangle, independent spacings | `2(nx+ny−2)` |

An unknown shape symbol raises an `ArgumentError` listing the valid options.

## The pairwise-geometry helper

[`borefield_geometry`](@ref) computes the geometry that spatial superposition operates on. Given the
coordinates and the borehole radius it returns the pairwise distance matrix `r` (with the diagonal
set to `rb`) and the flow-relative angle matrix `θ` in degrees (with the diagonal set to `0`). The
isotropic models (ILS, ICS, FLS) use only `r`; the moving models (MILS, MFLS) additionally use `θ`,
because under groundwater flow the pairwise response is asymmetric and depends on both the separation
and the angle to the flow. Because any two borehole pairs with the same `(r, θ)` give an identical
response, the distinct combinations of a layout, and how often each recurs, measure how much work
a solver can save by evaluating each unique geometry only once; the [`borefield_geometry`](@ref)
docstring shows how to compute those combinations directly from `r` and `θ`.

## Functions on this page

```@docs
borefield
borefield_geometry
borefield_rectangle
borefield_line
borefield_circle
borefield_L
borefield_U
borefield_open_rectangle
```
