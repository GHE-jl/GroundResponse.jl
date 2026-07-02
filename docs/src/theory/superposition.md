# Spatial superposition

A single-borehole g-function describes one isolated source. In a **borehole field** the boreholes
interact thermally: each one warms (or cools) the ground around its neighbours, so the field
response is larger than that of an isolated borehole. **Spatial superposition** combines the
single-borehole responses at every pairwise distance into one effective field g-function.

Both methods in this package assume the *Type II* boundary condition: the heat flux is uniform
along each borehole, and the **mean borehole-wall temperature is equal across all boreholes** (as
happens when the boreholes are connected in parallel to a common header). The output is again
normalised to a 1 W/m impulse.

## The pairwise-radius matrix

The raw ingredient is the ``n_b \times n_b`` matrix of distances between every pair of boreholes,
with the diagonal set to the borehole radius ``r_b`` (a borehole's distance "to itself" is its own
wall). [`borefield_radius`](@ref) builds this matrix from the coordinate array `xy`, along with the
flattened vector, the unique distances, an index map and the azimuth angles. Evaluating a model at
this matrix yields the ``n_t \times n_b \times n_b`` g-array that both superposition methods consume:

```math
g_{ij}(t) = g\bigl(r_{ij}, t\bigr), \qquad
r_{ij} = \begin{cases} \lVert \mathbf{x}_i - \mathbf{x}_j \rVert & i \ne j \\ r_b & i = j. \end{cases}
```

## Successive flux method

[`successive_flux`](@ref) (Nguyen & Pasquier, 2021) solves the field response **iteratively**. It
starts from a block-matrix estimate of the per-borehole fluxes that enforce equal mean temperatures,
then repeatedly:

1. forms the incremental (step) fluxes from the current flux estimate;
2. computes each borehole's temperature from the pairwise convolutions of those fluxes with the
   single-borehole responses ``g_{ij}`` (via FFT, [`DSP.jl`](https://github.com/JuliaDSP/DSP.jl));
3. measures the spread in borehole temperatures and corrects the fluxes to reduce it.

Iteration stops when the relative temperature spread is small or the residual stops improving. The
method is fast and memory-efficient and scales well from small to large fields, which makes it the
default used by [`ground_response`](@ref) for multi-borehole inputs. The initial linear solve falls
back to an SVD-based minimum-norm solution if the system is singular.

## Block matrix method

[`bloc_matrix`](@ref) (Dusseault et al., 2018) assembles the full space–time convolution system as a
single block matrix and solves it with **one direct linear solve**. Each ``(i, j)`` block is the
Toeplitz convolution operator built from ``g_{ij}``; an extra row/column enforces the equal-mean-
temperature constraint and a unit total load. The direct solve is more expensive in memory and time
than the successive-flux iteration, but it is unconditionally stable and serves as the reference
against which the iterative method is validated.

## Two ways to call each method

Both functions accept either a **precomputed g-array** (low-level) or **borefield parameters with a
model** (high-level):

```julia
# High-level: pass coordinates and a model — the radius matrix and g-array are built internally
g = successive_flux(t, rb, xy, FLSModel(150.0, 4.0, ks, Cs))
g = bloc_matrix(t, rb, xy, MILSModel(rb, ks, Cs, Cf, vD))

# Low-level: pass a precomputed nt × nb × nb g-array (any model, or a custom one)
r3D = borefield_radius(xy, rb)[1]          # nb×nb radius matrix
g3D = fls(t, r3D, 150.0, 4.0, ks, Cs)      # nt × nb × nb
g   = successive_flux(g3D)
```

The high-level overloads are exactly the extension point for custom models: a new
`AbstractGroundModel` becomes usable everywhere by defining its `successive_flux` /
`bloc_matrix` methods (see [Overview](@ref)).

## Functions on this page

```@docs
successive_flux
bloc_matrix
```
