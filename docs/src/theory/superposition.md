# Spatial superposition

A single-borehole g-function describes one isolated source. In a **borehole field** the boreholes
interact thermally: each one warms (or cools) the ground around its neighbours, so the field
response is larger than that of an isolated borehole. **Spatial superposition** combines the
single-borehole responses at every pairwise distance into one effective field g-function.

The output is always normalised to a unit total-field impulse and returned in °C·m/W, whatever the
boundary condition. As with the single-borehole models, the internal space–time solve produces the
borefield g-function `ḡ(t)` for a unit step in time; **load-side temporal superposition** (the
convolution of `ḡ` with an actual time-varying load) is applied separately, downstream of this
package.

## The three boundary conditions

How the heat load is allowed to distribute among — and along — the boreholes defines the *boundary
condition* of the superposition. The package provides the three classical choices, in order of
increasing detail:

| | Boundary condition | Physical assumption | Function(s) |
|---|---|---|---|
| **BC-I**  | Uniform heat flux | Every borehole emits the **same** flux; no coupling constraint | [`uniform_flux`](@ref) |
| **BC-II** | Uniform flux per borehole, equal mean temperature | Flux is uniform **along** each borehole; the **mean** wall temperature is equal across boreholes (parallel connection) | [`successive_flux`](@ref), [`bloc_matrix`](@ref) |
| **BC-III**| Uniform borehole-wall temperature | Flux varies **within** each borehole so the wall temperature is uniform over the **whole field** | [`segment_response`](@ref) |

All three are available and interchangeable through the same high-level interface; they bracket the
true field response, with `BC-I ≥ BC-II ≥ BC-III` at any given time (a finer boundary condition
lets the field shed heat more efficiently, lowering the g-function). BC-I and BC-II coincide for a
field whose boreholes are all geometrically equivalent (e.g. a symmetric 2×2 grid).

## The pairwise-geometry matrices

The raw ingredient of every method is the ``n_b \times n_b`` matrix of distances between every pair
of boreholes, with the diagonal set to the borehole radius ``r_b`` (a borehole's distance "to
itself" is its own wall). [`borefield_geometry`](@ref) builds this matrix from the coordinate array
`xy`, together with the flow-relative angle matrix ``\theta`` (in degrees, diagonal ``0``), the
unique ``(r, \theta)`` pairs, an index map and the per-pair counts. The isotropic models use only
``r``; the moving models (MILS, MFLS) also use ``\theta`` because advection makes the pairwise
response asymmetric. Evaluating a model at these matrices yields the ``n_t \times n_b \times n_b``
g-array that the superposition methods consume:

```math
g_{ij}(t) = g\bigl(r_{ij}, \theta_{ij}, t\bigr), \qquad
r_{ij} = \begin{cases} \lVert \mathbf{x}_i - \mathbf{x}_j \rVert & i \ne j \\ r_b & i = j, \end{cases}
```

where ``\theta_{ij}`` is the angle of the source→receiver direction relative to the flow (``+x``);
for the isotropic models the ``\theta`` dependence drops out.

## BC-I — uniform flux

[`uniform_flux`](@ref) (the M1 model of Guo et al., 2021) is the simplest scheme: every borehole
emits the same constant flux, and the field g-function is the average over all boreholes of each
borehole's self-plus-mutual response. It requires **no linear solve** — just a double sum — which
makes it numerically robust even for the strongly asymmetric response matrices produced by the
advection models. It is the upper bound of the three and the scheme Guo et al. (2021) validated
against a 3-D finite-element model with groundwater flow.

## BC-II — uniform flux per borehole, equal mean temperature

BC-II adds the coupling constraint that all boreholes reach the **same mean wall temperature** (as
happens when they are connected in parallel to a common header), while the flux is still uniform
*along* each borehole. Two equivalent solvers are provided:

[`successive_flux`](@ref) (Nguyen & Pasquier, 2021) solves the field response **iteratively**. It
starts from a block-matrix estimate of the per-borehole fluxes that enforce equal mean temperatures,
then repeatedly forms the incremental fluxes, computes each borehole's temperature from the pairwise
convolutions of those fluxes with the single-borehole responses ``g_{ij}`` (via FFT,
[`DSP.jl`](https://github.com/JuliaDSP/DSP.jl)), and corrects the fluxes to reduce the temperature
spread. It is fast and memory-efficient, scales well from small to large fields, and is the default
used by [`ground_response`](@ref) for multi-borehole inputs. The initial linear solve falls back to
an SVD-based minimum-norm solution if the system is singular.

Because that convolution is performed by FFT, it assumes a **constant time step**; a non-uniform
(e.g. log-spaced) `t` would corrupt the temporal superposition. The `interp = true` default
therefore evaluates the solver on an internal constant-step grid — the geometric hour→decade
sub-sampling of Nguyen & Pasquier (2021) — and interpolates back to the requested `t` with a
monotone cubic (PCHIP). For this temporal solver `interp` is a correctness requirement on any
non-uniform `t`, and it also bounds the cost to the sub-sample nodes regardless of `length(t)`; pass
`interp = false` only when `t` is already uniform. The same option is offered for the instantaneous
methods ([`bloc_matrix`](@ref), [`uniform_flux`](@ref)), where it is purely a performance choice.

[`bloc_matrix`](@ref) (Dusseault et al., 2018) assembles the full space–time convolution system as a
single block matrix and solves it with **one direct linear solve**. Each ``(i, j)`` block is the
Toeplitz convolution operator built from ``g_{ij}``; an extra row/column enforces the equal-mean-
temperature constraint and a unit total load. The direct solve is more expensive in memory and time
than the successive-flux iteration, but it is unconditionally stable and serves as the reference
against which the iterative method is validated.

## BC-III — uniform borehole-wall temperature (segments)

BC-III (Cimmino & Bernier, 2014) is the finest boundary condition. Each borehole is divided into
`nseg` **segments**, and the heat flux is allowed to vary from segment to segment so that the
borehole-wall temperature is uniform over the whole field and constant with depth. Building the
segment response requires the **segment-to-segment finite line source** — the general
[`fls`](@ref) kernel with independent emitting ``(H_1, D_1)`` and receiving ``(H_2, D_2)`` segments
— which reduces exactly to the classical FLS when the two coincide.

[`segment_response`](@ref) then **generalises the block-matrix formulation of BC-II from ``n_b``
boreholes to ``n_b \cdot n_{seg}`` segments**: it assembles one space–time block-Toeplitz system
whose unknowns are the segment fluxes and the single common wall temperature, with an energy-balance
row weighted by segment length. Because the flux profile along each borehole redistributes over
time to keep the wall isothermal, this internal solve intrinsically couples space and time — exactly
as `bloc_matrix` already does for BC-II. With `nseg = 1` it reduces identically to
[`bloc_matrix`](@ref). This coupling is internal to constructing `ḡ(t)` and is distinct from the
load-side temporal superposition applied downstream.

A second BC-III solver, [`segment_response_marching`](@ref), is also available. It marches forward
one time step at a time (a small per-step solve) and carries the load history by **incremental
temporal superposition**, following Cimmino (2018). It is far cheaper than the block matrix
(`O(nS·nt)` memory versus `O((nS·nt)²)`) and is the temporally rigorous formulation. Note that
`segment_response` (block matrix) and `segment_response_marching` are *different temporal
formulations*: the block matrix applies spatial superposition independently at each time step (no
load history), whereas the marching method superposes the history — the two agree only when the flux
never redistributes in time (`nseg = 1` on a symmetric field) and otherwise differ by a few percent.
See `dev/dev_bc3_marching_vs_block.jl` for the comparison.

!!! note "Scope and cost"
    BC-III is implemented for the finite line source (vertical boreholes). The direct block solve
    scales as ``(n_b \cdot n_{seg} \cdot n_t)`` unknowns, so memory and time grow quickly with
    field size and segment count; a handful of segments (`nseg = 6`–`12`) is usually enough for
    converged g-functions.

## Selecting the boundary condition

The boundary condition is chosen on the [`FLSModel`](@ref) through its optional `nseg` field, so the
same high-level call serves all three:

```julia
xy = borefield(:rectangle, 5, 5, 5.0)

g_bc1 = uniform_flux(t, rb, xy, FLSModel(150.0, 4.0, ks, Cs))        # BC-I
g_bc2 = ground_response(t, rb, xy, FLSModel(150.0, 4.0, ks, Cs))     # BC-II (default, nseg = 1)
g_bc3 = ground_response(t, rb, xy, FLSModel(150.0, 4.0, ks, Cs, 8))  # BC-III, 8 segments
```

## Two ways to call each method

Every method accepts either a **precomputed g-array** (low-level) or **borefield parameters with a
model** (high-level):

```julia
# High-level: pass coordinates and a model — the geometry and g-array are built internally
g = successive_flux(t, rb, xy, FLSModel(150.0, 4.0, ks, Cs))
g = bloc_matrix(t, rb, xy, MILSModel(rb, ks, Cs, Cf, vD))
g = segment_response(t, rb, xy, FLSModel(150.0, 4.0, ks, Cs, 8))

# Low-level: pass a precomputed g-array (any model, or a custom one)
r3D = borefield_geometry(xy, rb)[1]        # nb×nb distance matrix
g3D = fls(t, r3D, 150.0, 4.0, ks, Cs)      # nt × nb × nb
g   = successive_flux(g3D)
```

A single method is the extension point for custom models: a new `AbstractGroundModel` becomes usable
across every method — `uniform_flux`, `successive_flux`, `bloc_matrix`, and `ground_response` — by
defining `GroundResponse._response_array(t, rb, xy, m)`, which returns the pairwise `nt × nb × nb`
response array (see [Overview](@ref)).

## Functions on this page

```@docs
uniform_flux
successive_flux
bloc_matrix
segment_response
segment_response_marching
```
