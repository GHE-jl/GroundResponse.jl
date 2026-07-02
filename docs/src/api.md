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

- [`successive_flux`](@ref) — iterative method (Nguyen & Pasquier, 2021)
- [`bloc_matrix`](@ref) — direct block-matrix solve (Dusseault et al., 2018)

### Borefield layouts

- [`borefield`](@ref) — unified entry point
- [`borefield_radius`](@ref) — pairwise-distance matrix
- [`borefield_rectangle`](@ref), [`borefield_line`](@ref), [`borefield_circle`](@ref)
- [`borefield_L`](@ref), [`borefield_U`](@ref), [`borefield_open_rectangle`](@ref)

!!! tip "Where the docstrings live"
    Full signatures and argument lists are rendered inline on the theory pages:
    [Overview](@ref), [Line-source models](@ref), [Moving-source models](@ref),
    [Spatial superposition](@ref) and [Borefields](@ref).

!!! note "Short-term ANN model"
    The package also contains a short-term g-function based on an artificial neural network
    (`gST_ANN`, Pasquier et al. 2018). It is not yet fully integrated and is therefore **not
    exported**; it does not appear in the index above.
