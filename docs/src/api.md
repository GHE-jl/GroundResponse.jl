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

### Short-term ANN model

- [`short_term_response`](@ref) — ANN-based short-term outlet transfer function
  (Pasquier, Zarrella & Labib, 2018)
- [`short_term_nodes`](@ref) — raw ANN transfer function on its 85 native time nodes (the
  building block behind `short_term_response` and `outlet_transfer_function`)

### Spatial superposition methods

Presented in order of increasing boundary-condition detail (BC-I → BC-II → BC-III):

- [`uniform_flux`](@ref) — **BC-I**: equal flux on every borehole (Guo et al., 2021)
- [`successive_flux`](@ref) — **BC-II**: iterative equal-mean-temperature solve (Nguyen & Pasquier, 2021)
- [`bloc_matrix`](@ref) — **BC-II**: direct block-matrix solve (Dusseault et al., 2018)
- [`segment_response`](@ref) — **BC-III**: segment (uniform-wall-temperature) solve (Cimmino & Bernier, 2014)
- [`segment_response_marching`](@ref) — **BC-III**: stepwise time-marching alternative to the block
  matrix (Cimmino, 2018)

The boundary condition is selected on [`FLSModel`](@ref) through its optional `nseg` field
(`nseg > 1` ⇒ BC-III).

### Borefield layouts

- [`borefield`](@ref) — unified entry point
- [`borefield_geometry`](@ref) — pairwise distance and angle matrices
- [`borefield_rectangle`](@ref), [`borefield_line`](@ref), [`borefield_circle`](@ref)
- [`borefield_L`](@ref), [`borefield_U`](@ref), [`borefield_open_rectangle`](@ref)

!!! tip "Where the docstrings live"
    Full signatures and argument lists are rendered inline on the theory pages:
    [Overview](@ref), [Line-source models](@ref), [Moving-source models](@ref),
    [Spatial superposition](@ref) and [Borefields](@ref). The short-term ANN model is documented
    below since it does not yet have a dedicated theory page.

## Short-term ANN model

```@docs
short_term_response
short_term_nodes
```
