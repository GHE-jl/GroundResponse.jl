# GroundResponse.jl — Claude Code Context

## Project overview

Julia package implementing ground thermal response models for borehole heat exchanger (BHE)
design and analysis. Models: ILS, ICS, FLS (static), MILS, MFLS (with groundwater advection).

Key directories:
- `src/` — model implementations + spatial superposition logic
- `script/` — demonstration and comparison scripts
- `docs/` — (new, not yet populated)
- `.github/` — CI config (new)

## Active investigation — MILS spatial superposition bug

### Symptom
Running `successive_flux(t, rb, xy33, MILSModel(...))` on a 3×3 borefield produces
wildly diverging g-function values (~±4000 °Cm/W) instead of values near zero like ILS/FLS.
See `script/script_spatial_superposition.jl:54` (has a `#TODO: Fix for the MILS` comment).

### Root cause

The article (Pasquier & Lamarche 2022, `src/moving_infinite_line_source.jl`) only provides
the **azimuthal mean** temperature `ḡ(r, t)` (Eq. 3) — a function of scalar distance `r` only.
This is what `_mils_series` computes.

The code extends this to a **directional** formula (Eq. 1):

    θ(x, y, t) = ḡ(r, t) × exp(x·vT/2α) / I₀(r·vT/2α)

This is mathematically correct and is implemented in `_mils` (used for 2D spatial maps in
`script_groundwater_advection.jl`).

**The problem**: when `mils(t, xy_matrix, rb, ...)` is called for spatial superposition, it
builds the 3D g-matrix using `_mils(dx, dy)` where `(dx, dy) = xy[i,:] - xy[j,:]`. Because
the directional formula depends on the sign of `dx`, the resulting matrix is **asymmetric**:
`g[i,j] ≠ g[j,i]`.

With `vD = 1e-6 m/s`, `Cf = 4.18e6`, `ks = 3.0`, `Cs = 2e6`:

    vT/(2α) = vD·Cf / (2·ks) ≈ 0.697 m⁻¹

For a borehole pair with `dx = ±10 m` (opposite corners of a 3×3 borefield, B=5 m):

    g[downstream] / g[upstream] = exp(2 × 10 × 0.697) ≈ 1.1 × 10⁶

The `successive_flux` and `bloc_matrix` algorithms assume a **symmetric** g-matrix (which
holds for ILS/ICS/FLS since response depends only on `r`). The 10⁶ asymmetry causes
the iterative solver to diverge.

### Proposed fix

In `src/spatial_superposition.jl`, change the `MILSModel` overloads of `successive_flux`
and `bloc_matrix` to use `_mils_series` (azimuthal mean) evaluated at the symmetric
inter-borehole distances — exactly as ILS does via `borefield_radius`:

```julia
function successive_flux(t, rb, xy, m::MILSModel)
    @assert all(>(0), t) "t must be strictly positive"
    r, _, _, _, _, _ = borefield_radius(xy, rb)   # symmetric distance matrix
    nb = size(r, 1); nt = length(t)
    T = float(promote_type(eltype(t), typeof(m.ks), typeof(m.Cs), typeof(m.Cf), typeof(m.vD)))
    g3D = Array{T,3}(undef, nt, nb, nb)
    for j in 1:nb, i in 1:nb
        ri = T(r[i, j])
        for k in eachindex(t)
            g3D[k, i, j] = _mils_series(T(t[k]), T(m.ks), T(m.Cs), T(m.Cf), ri, T(m.vD))
        end
    end
    return successive_flux(g3D)
end
```

Apply the same change to `bloc_matrix(t, rb, xy, m::MILSModel)`.

### Why azimuthal mean is correct for spatial superposition

The Type II boundary condition (all boreholes at same mean temperature, used by both
`successive_flux` and `bloc_matrix`) requires averaging the temperature over the borehole
wall of borehole `i` due to source `j`. Since borehole `i`'s wall is a small circle at
distance `r_ij` from source `j`, this average is `_mils_series(t, r_ij, ...)`.
The directional `_mils` formula is correct for **point evaluation** (2D maps, TRT analysis)
but not for wall-averaged borefield g-functions.

### Files involved
- `src/moving_infinite_line_source.jl` — `_mils_series` (azimuthal mean), `_mils` (directional)
- `src/spatial_superposition.jl:103-109` — `successive_flux` and `bloc_matrix` MILS overloads
- `script/script_spatial_superposition.jl:54` — triggers the bug, has TODO comment

### Reference
Pasquier, P., & Lamarche, L. (2022). Analytic expressions for the moving infinite line source
model. Geothermics, 103, 102413. https://doi.org/10.1016/j.geothermics.2022.102413
PDF: `C:\Users\gabri\OneDrive\Recherche\Zotero\storage\7UH3Y79T\`
