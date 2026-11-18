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

### Borefield response arrays

- [`gmatrix`](@ref) — pairwise response array of a borefield, per borehole or per segment

The spatial superposition of these arrays under BC-I to BC-IV (`bc_I`, `bc_II`, `bc_III`, `bc_IV`
and `ground_response`) is in
[GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl).

### Borefield layouts

- [`borefield`](@ref) — unified entry point
- [`borefield_geometry`](@ref) — pairwise distance and angle matrices
- [`borefield_rectangle`](@ref), [`borefield_line`](@ref), [`borefield_circle`](@ref)
- [`borefield_L`](@ref), [`borefield_U`](@ref), [`borefield_open_rectangle`](@ref)

!!! tip "Where the docstrings live"
    Full signatures and argument lists are rendered inline on the theory pages:
    [Overview](@ref), [Line-source models](@ref), [Moving-source models](@ref) and
    [Borefields](@ref).
