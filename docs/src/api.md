# API reference

This page is an alphabetical index of every exported symbol. Each entry links to its full
docstring, which lives next to the relevant theory on the modeling pages.

```@index
Modules = [GroundResponse]
```

## By topic

### Ground model types

- [`AbstractGroundModel`](@ref) — the abstract supertype and extension point
- [`ILSModel`](@ref), [`ICSModel`](@ref), [`FLSModel`](@ref) — conductive models
- [`MILSModel`](@ref), [`MFLSModel`](@ref) — moving (advection) models

### Single-borehole g-functions

- [`ils`](@ref) — infinite line source
- [`ics`](@ref) — infinite cylindrical source
- [`fls`](@ref) — finite line source
- [`mils`](@ref) — moving infinite line source
- [`mfls`](@ref) — moving finite line source

### High-level interface

- [`ground_response`](@ref) — dispatch over model type and field size

### Spatial superposition methods

One function per boundary condition, in order of increasing detail (BC-I → BC-II → BC-III):

- [`bc_I`](@ref) — **BC-I**: equal flux on every borehole (Guo et al., 2021)
- [`bc_II`](@ref) — **BC-II**: uniform flux per borehole, equal mean wall temperature, solved by
  time marching
- [`bc_III`](@ref) — **BC-III**: segments with a uniform wall temperature (Cimmino & Bernier, 2014),
  the same marching solver

[`ground_response`](@ref) selects them with its `bc` keyword; BC-III is the default for an
[`FLSModel`](@ref) with `nseg > 1`.

Legacy solvers, superseded by `bc_II` / `bc_III` and kept for comparison (the names without
`_legacy` are deprecated aliases):

- [`successive_flux_legacy`](@ref) — BC-II, iterative (Nguyen & Pasquier, 2021)
- [`bloc_matrix_legacy`](@ref) — BC-II, block matrix (Dusseault et al., 2018)
- [`segment_response_legacy`](@ref) — BC-III, block matrix
- [`segment_response_marching_legacy`](@ref) — BC-III, time marching on strided slices

### Borefield layouts

- [`borefield`](@ref) — unified entry point
- [`borefield_geometry`](@ref) — pairwise distance and angle matrices
- [`borefield_rectangle`](@ref), [`borefield_line`](@ref), [`borefield_circle`](@ref)
- [`borefield_L`](@ref), [`borefield_U`](@ref), [`borefield_open_rectangle`](@ref)

!!! tip "Where the docstrings live"
    Full signatures and argument lists are rendered inline on the theory pages:
    [Overview](@ref), [Line-source models](@ref), [Moving-source models](@ref),
    [Spatial superposition](@ref) and [Borefields](@ref).
