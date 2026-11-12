# Spatial superposition

A single-borehole g-function describes one isolated source. In a **borehole field** the boreholes
interact thermally: each one warms (or cools) the ground around its neighbours, so the field
response is larger than that of an isolated borehole. **Spatial superposition** combines the
single-borehole responses at every pairwise distance into one effective field g-function.

The output is always normalised to 1 W/m on every borehole and returned in °C·m/W, whatever the
boundary condition. The internal space–time solve produces the borefield g-function `ḡ(t)` for a
unit step in time. **Load-side temporal superposition** (the convolution of `ḡ` with an actual
time-varying load) is applied separately, downstream of this package.

## The three boundary conditions

How the heat load is allowed to distribute among, and along, the boreholes defines the *boundary
condition* of the superposition. The package provides the three classical choices, in order of
increasing detail, each with its own function:

| | Boundary condition | Physical assumption | Function |
|---|---|---|---|
| **BC-I**  | Uniform heat flux | Every borehole emits the **same** flux; no coupling constraint | [`bc_I`](@ref) |
| **BC-II** | Uniform flux per borehole, equal mean temperature | Flux is uniform **along** each borehole; the **mean** wall temperature is equal across boreholes (parallel connection) | [`bc_II`](@ref) |
| **BC-III**| Uniform borehole-wall temperature | Flux varies **within** each borehole so the wall temperature is uniform over the **whole field** | [`bc_III`](@ref) |

They bracket the true field response, with `BC-I ≥ BC-II ≥ BC-III` at any given time (a finer
boundary condition lets the field shed heat more efficiently, lowering the g-function). BC-I and
BC-II coincide for a field whose boreholes are all geometrically equivalent (e.g. a symmetric 2×2
grid), since the fluxes then stay equal.

Under BC-II and BC-III, the 1 W/m is the **average** over the field: the flux of each borehole (or
segment) is an unknown, different from one borehole to the next and changing in time. A corner
borehole, for instance, takes a larger share of the load than a central one once the boreholes
interact. That is why these two conditions superpose the load history in time as well as in space.

## The pairwise-geometry matrices

The raw ingredient of every method is the ``n_b \times n_b`` matrix of distances between every pair
of boreholes, with the diagonal set to the borehole radius ``r_b`` (a borehole's distance "to
itself" is its own wall). [`borefield_geometry`](@ref) builds this matrix from the coordinate array
`xy`, together with the flow-relative angle matrix ``\theta`` (in degrees, diagonal ``0``). The
isotropic models use only ``r``, and the moving models (MILS, MFLS) also use ``\theta`` because
advection makes the pairwise response asymmetric. Evaluating a model at these matrices yields the
``n_t \times n_b \times n_b`` g-array that the superposition methods consume:

```math
g_{ij}(t) = g\bigl(r_{ij}, \theta_{ij}, t\bigr), \qquad
r_{ij} = \begin{cases} \lVert \mathbf{x}_i - \mathbf{x}_j \rVert & i \ne j \\ r_b & i = j, \end{cases}
```

where ``\theta_{ij}`` is the angle of the source→receiver direction relative to the flow (``+x``);
for the isotropic models the ``\theta`` dependence drops out.

## BC-I — uniform flux

[`bc_I`](@ref) (the M1 model of Guo et al., 2021) is the simplest scheme: every borehole emits the
same constant flux, and the field g-function is the average over all boreholes of each borehole's
self-plus-mutual response:

```math
\bar g(t) = \frac{1}{n_b} \sum_{i=1}^{n_b} \sum_{j=1}^{n_b} g_{ij}(t).
```

It requires **no linear solve**, which makes it numerically robust even for the strongly asymmetric
response matrices produced by the advection models. It is the upper bound of the three and the
scheme Guo et al. (2021) validated against a 3-D finite-element model with groundwater flow.

## BC-II and BC-III — time marching

BC-II adds the coupling constraint that all boreholes reach the **same mean wall temperature** (as
happens when they are connected in parallel to a common header), while the flux is still uniform
*along* each borehole. BC-III (Cimmino & Bernier, 2014) divides each borehole into `nseg`
**segments** and lets the flux vary from segment to segment, so that the wall temperature is uniform
over the whole field and constant with depth. With `nseg = 1` the two are the same condition, and
they are solved by the same algorithm over ``n_S`` sources: the ``n_b`` boreholes for
[`bc_II`](@ref), or the ``n_b \cdot n_{seg}`` segments for [`bc_III`](@ref). The segment responses
use the **segment-to-segment finite line source**, the general [`fls`](@ref) kernel with independent
emitting ``(H_1, D_1)`` and receiving ``(H_2, D_2)`` segments.

The flux of each source is constant over each time step, so the unknowns of step ``p`` are the flux
increments ``f_{j,p} = q_j(t_p) - q_j(t_{p-1})`` and the common wall temperature ``T_p``. The
temperature of source ``i`` superposes the response to every increment since ``t = 0``:

```math
T_i(t_p) = \sum_{j} \sum_{l=1}^{p} g_{ij}(t_p - t_{l-1})\, f_{j,l}
         = \sum_{j} g_{ij}(t_p - t_{p-1})\, f_{j,p} + h_{i,p},
\qquad
h_{i,p} = \sum_{l<p} \sum_{j} g_{ij}(t_p - t_{l-1})\, f_{j,l},
```

where the load history ``h_{i,p}`` only involves the increments of earlier steps, already known.
With one common temperature (``T_i(t_p) = T_p`` for every ``i``) and a constant total heat rate
(``\sum_j H_j f_{j,p} = H_{ref}`` at ``p = 1``, ``0`` afterwards), each step solves a small
``(n_S + 1) \times (n_S + 1)`` system:

```math
\begin{bmatrix} g(t_p - t_{p-1}) & -\mathbf{1} \\ \mathbf{H}^\top & 0 \end{bmatrix}
\begin{bmatrix} \mathbf{f}_p \\ T_p \end{bmatrix}
=
\begin{bmatrix} -\mathbf{h}_p \\ H_{ref}\,\delta_{p1} \end{bmatrix}.
```

The g-function is ``\bar g(t_p) = T_p`` (with ``H_{ref}`` the total length, i.e. 1 W/m per
borehole). This is a **direct** solution: no iteration and no tolerance. Stacking every step gives the
single block system of Dusseault et al. (2018), which is lower block-triangular in time; time
marching is its forward substitution, so both give the same solution, at a cost of
``n_t\, n_S^3`` (one factorisation) plus ``n_t^2 n_S^2 / 2`` for the history instead of
``(n_t\, n_S)^3``. On a uniform grid ``t_k = k\,\Delta t`` the lags are themselves grid times,
``g(t_p - t_{l-1}) = g(t_{p-l+1})``, so only ``g`` at the requested times is needed and the step
matrix ``g(\Delta t)`` is factorised once.

The forward substitution multiplies rounding errors by ``g(t_k)/g(\Delta t)`` at every step. It is
therefore unstable when the first step lies in the onset of the self response
(``\Delta t \lesssim r_b^2/\alpha``, below about one hour for usual boreholes), where ``g`` grows by
orders of magnitude from one step to the next. Heat has then not left the vicinity of the borehole
wall: the sources do not interact and their fluxes stay uniform, so the uniform-flux (BC-I) value is
returned for such a block (detected when a self response grows more than tenfold over the first
step). The legacy BC-III solvers have no such safeguard and are wrong on those blocks.

## Time sampling

Both BC-II and BC-III therefore need a uniform grid ``t = \Delta t, 2\Delta t, \dots, n\Delta t``.
With `interp = true` (the default), the g-function is solved on internal uniform blocks — hourly to
one day, daily to one month, monthly to one year, quarterly to ten years and yearly to fifty years
(decades beyond) — following the geometric progression of Cimmino et al. (2013) and Nguyen & Pasquier
(2021), and the nodes are interpolated to the requested `t` by a monotone cubic (PCHIP) in
``\log t``. This keeps the cost to about two hundred solved times whatever the length of `t`. Each
block restarts from ``t = 0`` with its own step, so the flux redistribution during its first step is
not resolved; the quarterly and yearly blocks keep that error to about 0.2 % for a 100-borehole
field. Pass `interp = false` to solve directly on a `t` that is already ``\Delta t, 2\Delta t, \dots``
(any other `t` is rejected). For BC-I `interp` is a performance choice only.

!!! note "Scope and cost"
    BC-III is implemented for the finite line source (vertical boreholes). Memory grows as
    ``n_t (n_b\, n_{seg})^2``, so it grows quickly with field size and segment count; a handful of
    segments (`nseg = 2`–`12`) is usually enough for converged g-functions.

## Selecting the boundary condition

[`ground_response`](@ref) selects the boundary condition with its `bc` keyword; the number of
segments is a field of [`FLSModel`](@ref), and BC-III is the default when it is larger than one:

```julia
xy = borefield(:rectangle, 5, 5, 5.0)

g_bc1 = ground_response(t, rb, xy, FLSModel(150.0, 4.0, ks, Cs); bc = :I)   # BC-I
g_bc2 = ground_response(t, rb, xy, FLSModel(150.0, 4.0, ks, Cs))            # BC-II (nseg = 1)
g_bc3 = ground_response(t, rb, xy, FLSModel(150.0, 4.0, ks, Cs, 8))         # BC-III, 8 segments
```

## Two ways to call each method

Every method accepts either a **precomputed g-array** (low-level) or **borefield parameters with a
model** (high-level):

```julia
# High-level: pass coordinates and a model — the geometry and g-array are built internally
g = bc_I(t, rb, xy, MILSModel(rb, ks, Cs, Cf, vD))
g = bc_II(t, rb, xy, FLSModel(150.0, 4.0, ks, Cs))
g = bc_III(t, rb, xy, FLSModel(150.0, 4.0, ks, Cs, 8))

# Low-level: pass a precomputed g-array (any model, or a custom one) on t = Δt, 2Δt, …
tu  = 2_592_000.0 .* (1:12)                # monthly steps
r   = borefield_geometry(xy, rb)[1]        # nb × nb distance matrix
g3D = fls(tu, r, 150.0, 4.0, ks, Cs)       # nt × nb × nb
g   = bc_II(g3D)
```

A single method is the extension point for custom models: a new `AbstractGroundModel` becomes usable
across `bc_I`, `bc_II` and `ground_response` by defining
`GroundResponse._borehole_response(t, rb, xy, m)`, which returns the pairwise `nt × nb × nb` response
array (see [Overview](@ref)).

## Legacy solvers

The solvers used before version 2.2 are kept under a `_legacy` suffix for comparison; their former
names still work as deprecated aliases. All of them need the same uniform time grid (provided by
`interp = true`).

- [`successive_flux_legacy`](@ref) (BC-II, Nguyen & Pasquier, 2021) iterates on the fluxes of every
  borehole at every time until the wall temperatures agree within 0.1 %, with one convolution per
  pair of boreholes. It converges towards the [`bc_II`](@ref) solution, within about ``10^{-5}`` to
  ``3 \times 10^{-4}``.
- [`bloc_matrix_legacy`](@ref) (BC-II, Dusseault et al., 2018) and
  [`segment_response_legacy`](@ref) (BC-III) solve the whole space–time block system at once. They
  give the same solution as [`bc_II`](@ref) and [`bc_III`](@ref), at a cost growing as
  ``(n_t\, n_S)^3``.
- [`segment_response_marching_legacy`](@ref) (BC-III) is the marching algorithm reading each lag of
  the load history as a strided slice of `g`, which is slower on large problems.

!!! note "Correction in 2.2"
    Before version 2.2, the block solvers (`bloc_matrix`, `segment_response`) filled column ``l`` of
    each time block with ``g(t_k)`` instead of ``g(t_k - t_{l-1})``, which drops the load history:
    each time step was solved as if the current fluxes had been applied since ``t = 0``. The error
    was zero on a symmetric field and reached a few percent once the fluxes redistribute. The legacy
    block solvers are corrected.

## Functions on this page

```@docs
bc_I
bc_II
bc_III
successive_flux_legacy
bloc_matrix_legacy
segment_response_legacy
segment_response_marching_legacy
```
