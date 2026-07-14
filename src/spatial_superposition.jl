using LinearAlgebra
using DSP: conv

# -------------------------------------------------------------------------------------------------
# Model response array and time sub-sampling
# -------------------------------------------------------------------------------------------------

# Pairwise (nt × nb × nb) response array of a borefield under any ground model. Centralises the
# per-model evaluation so the spatial-superposition backends stay model-agnostic. The isotropic
# models (ILS/ICS/FLS) key their geometry on the passed borehole radius `rb`; the moving models
# (MILS/MFLS) carry their own `m.rb` and additionally need the flow-relative angle matrix.
_response_array(t, rb, xy, m::ILSModel) = ils(t, borefield_geometry(xy, rb)[1], m.ks, m.Cs)
_response_array(t, rb, xy, m::ICSModel) = ics(t, borefield_geometry(xy, rb)[1], rb, m.ks, m.Cs)
_response_array(t, rb, xy, m::FLSModel) = fls(t, borefield_geometry(xy, rb)[1], m.H, m.D, m.ks, m.Cs)
function _response_array(t, rb, xy, m::MILSModel)
    r, θ = borefield_geometry(xy, m.rb)
    return mils(t, r, θ, m.rb, m.ks, m.Cs, m.Cf, m.vD)
end
function _response_array(t, rb, xy, m::MFLSModel)
    r, θ = borefield_geometry(xy, m.rb)
    return mfls(t, r, θ, m.H, m.rb, m.D, m.ks, m.Cs, m.Cf, m.vD)
end

# Fixed sub-sampling scheme for the TEMPORAL spatial-superposition solvers (`successive_flux`,
# `segment_response_marching`). Their spectral / index-based convolution assumes a CONSTANT time
# step, so passing an arbitrary (e.g. log-spaced) `t` silently corrupts the temporal superposition.
# The g-function is instead solved on the blocks below — each UNIFORMLY spaced, so the convolution
# within a block is valid — following the geometric hour→decade progression of Nguyen & Pasquier
# (2021), with two sub-minute blocks prepended so a short-term model can resolve the first minute.
# Each block is solved independently, the results are stitched (a finer block wins at a shared
# boundary), and PCHIP-interpolated to the requested `t`. Blocks entirely outside the requested
# [t0, tend] window are skipped, so the cost is bounded regardless of `length(t)`. A caller who wants
# the response on their own grid can bypass this with `interp = false` and accept the constant-step
# requirement, or call the low-level array kernel directly. The same routine also serves the
# instantaneous / direct paths (`uniform_flux`, `bloc_matrix`, single-borehole models), where it is a
# pure performance approximation rather than a correctness requirement.
const _SUBSAMPLE_BLOCKS = ((2.0,          30),   # every 2 s   → 60 s
                           (120.0,        30),   # every 2 min → 1 h
                           (3600.0,       24),   # every 1 h   → 1 day
                           (86_400.0,     30),   # every 1 day → 1 month (30-day)
                           (2_592_000.0,  12),   # every 1 month → 1 year (360-day)
                           (31_104_000.0, 10))   # every 1 year → 1 decade
const _SUBSAMPLE_DECADE = 311_040_000.0          # 1 decade [s]; count grows to cover the request

function _append_strictly_increasing!(ts, gs, tb, gb)
    @inbounds for k in eachindex(tb)
        if isempty(ts) || tb[k] > ts[end] * (1 + 1e-12)
            push!(ts, tb[k]); push!(gs, gb[k])
        end
    end
    return nothing
end

# Solve a temporal g-function on the fixed sub-sampled grid, then PCHIP-interpolate to `t`.
# `solve(tblock)` must return the g-function evaluated on the uniform vector `tblock`.
function _subsampled_gfunction(t, solve)
    t0, tend = float(minimum(t)), float(maximum(t))
    ts = Float64[]; gs = Float64[]
    # Finest block to include: the coarsest whose step ≤ t0. This adapts the lower bound to the
    # request — it yields one node at/below t0 for a safe interpolation without resolving finer than
    # asked. Over-resolving below t0 is not just wasteful: for the FLS family below its characteristic
    # time the sub-minute blocks make the marching solve ill-posed (near-singular per-step matrix).
    steps  = ntuple(i -> _SUBSAMPLE_BLOCKS[i][1], length(_SUBSAMPLE_BLOCKS))
    finest = t0 < steps[1] ? steps[1] : maximum(s for s in steps if s ≤ t0)
    for (step, n) in _SUBSAMPLE_BLOCKS
        step < finest && continue         # finer than the requested t0 needs
        step > tend  && break             # step-ordered: block starts after the window
        tb = collect(range(step, step * n; length = n))
        _append_strictly_increasing!(ts, gs, tb, solve(tb))
    end
    if tend > _SUBSAMPLE_DECADE            # beyond 1 decade → add decade-spaced points
        ndec = max(2, ceil(Int, tend / _SUBSAMPLE_DECADE))
        tb = collect(range(_SUBSAMPLE_DECADE, ndec * _SUBSAMPLE_DECADE; length = ndec))
        _append_strictly_increasing!(ts, gs, tb, solve(tb))
    end
    @assert length(ts) ≥ 2 "sub-sampling produced too few nodes; pass interp = false"
    itp = Interpolator(ts, gs)
    return [itp(clamp(tt, ts[1], ts[end])) for tt in t]
end

# Optionally route a backend `solve(t)` through the sub-sampling grid. `interp = true` computes on the
# reduced constant-step grid and interpolates to `t`; `interp = false` computes directly on `t`. Used
# by `ground_response` to give one uniform keyword across every backend.
_apply_interp(t, interp::Bool, solve) = interp ? _subsampled_gfunction(t, solve) : solve(t)

"""
    uniform_flux(g)
    uniform_flux(t, rb, xy, m::AbstractGroundModel)

Spatial superposition under a **uniform, equal heat flux on every borehole** — the M1 model of
Guo et al. (2021), equivalent to boundary condition Type-I. 
Every borehole emits the same constant flux; the response of borehole `i` is the sum of its
self response and the mutual responses from all other boreholes, and the borefield g-function is
the average of these responses over all boreholes (Guo et al. 2021, Eqs. 11–12):

    Θᵢ = Σⱼ g[i, j]              (self j = i plus mutual j ≠ i)
    ḡ = (1 / nb²)Σᵢ Θᵢ

Contrary to `successive_flux` / `bloc_matrix` — which enforce a uniform mean borehole-wall
temperature (Type-II) through a linear solve — M1 requires no matrix inversion. It is therefore
numerically robust even for the strongly asymmetric, direction-dependent response matrices produced
by the advection models (`MILSModel`, `MFLSModel`). This is the scheme Guo et al. (2021) validated
against a 3D finite-element model for borefields with groundwater flow.

# Arguments
    - `g`: 3D response matrix (nt × nb × nb) or 2D matrix (nb × nb) for a single time. `g[k, i, j]`
        is the response at borehole `i` to a unit flux at borehole `j`. May be asymmetric.
    - `t`: Time vector (nt x 1) [s]
    - `rb`: Borehole radius [m]
    - `xy`: Matrix of borehole coordinates where the line source is at (0,0) (nb x 2) [m]
        - Can be computed with `borefield(:rectangle, ...)` from borefield.jl.
    - `m`: Ground model parameters (e.g. `FLSModel(150, 4, 3.0, 2e6)`)
# Output
    - `ḡ`: g-function of the borefield spatial superposition [°C·m/W]
# Reference
    - Guo, Y., Hu, X., Banks, J., & Liu, W. V. (2021). Considering buried depth for vertical
        borehole heat exchangers in a borehole field with groundwater flow — An extended solution.
        Energy and Buildings, 235, 110722. https://doi.org/10.1016/j.enbuild.2021.110722
"""
function uniform_flux(g::AbstractMatrix{<:Real})
    nb1, nb2 = size(g)
    @assert nb1 == nb2 "g must be nb × nb"
    @assert all(isfinite, g) "g must contain only finite values"
    return sum(g) / nb1^2
end
function uniform_flux(g::AbstractArray{<:Real,3})
    nt, nb1, nb2 = size(g)
    @assert nb1 == nb2 "g must be nt × nb × nb"
    @assert all(isfinite, g) "g must contain only finite values"
    return [sum(@view g[k, :, :]) / nb1^2 for k in 1:nt]
end
# Generic model-level method — routes through `_response_array`, so it serves the isotropic models
# (ILS/ICS/FLS) and any custom `AbstractGroundModel`. The moving models keep their own methods below
# to enforce t > 0 (their kernel is singular at zero time).
uniform_flux(t, rb, xy, m::AbstractGroundModel) = uniform_flux(_response_array(t, rb, xy, m))
function uniform_flux(t, rb, xy, m::MILSModel)
    @assert all(>(0), t) "t must be strictly positive to avoid a singular zero-time response"
    r, θ = borefield_geometry(xy, m.rb)
    return uniform_flux(mils(t, r, θ, m.rb, m.ks, m.Cs, m.Cf, m.vD))
end
function uniform_flux(t, rb, xy, m::MFLSModel)
    @assert all(>(0), t) "t must be strictly positive to avoid a singular zero-time response"
    r, θ = borefield_geometry(xy, m.rb)
    return uniform_flux(mfls(t, r, θ, m.H, m.rb, m.D, m.ks, m.Cs, m.Cf, m.vD))
end

"""
    successive_flux(g)
    successive_flux(t, rb, xy, m::AbstractGroundModel; interp = true)

Iteratively solve spatial superposition for a borefield using the successive flux approach of
Nguyen and Pasquier (2021) to obtain the g-functions of a borefield. This approach assumes that heat
flux is uniform along all the borehole, and that the mean temperature is the same for all
boreholes (Type II). The g-function generated is for an impulse of 1 W/m.

This is a **temporal** solver: it reproduces the load history by a spectral convolution, which
assumes a **constant time step**. The `interp` keyword (default `true`) solves on an internal
constant-step grid — the geometric hour→decade sub-sampling of Nguyen & Pasquier (2021) — and
PCHIP-interpolates to the requested `t`. For this solver `interp = true` is a **correctness**
requirement on any non-uniform (e.g. log-spaced) `t`, not merely a speed-up: it is the only way the
constant-step convolution stays valid. It also bounds the cost — the work is set by the ~one hundred
sub-sample nodes, independent of `length(t)`. Pass `interp = false` only when `t` is already
uniformly spaced, to solve directly on it (the caller then owns the constant-step requirement).
# Arguments
    - `g`: A 3D g-function matrix for all radius of the borefield (nt x nb x nb) [°Cm/W]
        - Each time step (1 x nb x nb) has the borefield response for each radius between boreholes.
        - Can be pre-computed with any ground model and passed directly.
    - `t`: Time vector (nt x 1) [s]
    - `rb`: Borehole radius [m]
    - `xy`: Matrix of borehole coordinates where the line source is at (0,0) (nb x 2) [m]
        - Can be computed with `borefield(:rectangle, ...)` from borefield.jl.
    - `m`: Ground model parameters (e.g. `FLSModel(150, 4, 3.0, 2e6)`)
# Output
    - `ḡ`: g-function of the borefield spatial superposition [°C·m/W]
# Reference
    - Nguyen, A., & Pasquier, P. (2021). A successive flux estimation method for rapid g-function
        construction of small to large-scale ground heat exchanger. Renewable Energy, 165, 359–368.
        https://doi.org/10.1016/j.renene.2020.10.074
"""
function successive_flux(g::AbstractArray{<:Real,3})
    # Basic parameters
    nt, nb1, nb2 = size(g)
    @assert nb1 == nb2 "g must be nt × nb × nb"
    @assert all(isfinite, g) "g must contain only finite values"
    nb = nb1

    # First estimation of g-function using block matrix (Eq. 20)
    GG = zeros(eltype(g), nt, nb + 1, nb + 1)
    @views GG[:, 1:nb, 1:nb] .= g
    @views GG[:, nb + 1, 1:nb] .= 1
    @views GG[:, 1:nb, nb + 1] .= 1

    b = zeros(eltype(g), nb + 1)
    b[end] = 1

    x = zeros(eltype(g), nt, nb)
    for it in 1:nt
        M = Matrix(@view GG[it, :, :])
        sol = try
            M \ b
        catch err
            if err isa SingularException
                @warn "Singular initial flux system at time index $it
                    Using an SVD-based minimum-norm solve instead."
                svd(M) \ b
            else
                rethrow()
            end
        end
        x[it, :] .= sol[1:nb]
    end

    # Successive flux estimation
    e1 = 10.0
    e2 = Inf
    e3 = Inf
    k  = 0
    kmax = 100
    ḡ = zeros(eltype(g), nt)

    while e3 > 0.15 && e1 > 1e-3 && e1 < e2 && k < kmax
        k += 1                                  # Iteration counter
        f = vcat(x[1:1, :], diff(x, dims=1))    # Step fluxes (Eq. 4)
        # Temperature response via pairwise convolutions (Eq. 6, 8)
        hh = zeros(eltype(g), nt, nb)
        for j in 1:nb, i in 1:nb
            hh[:, i] .+= conv(f[:, j], g[:, i, j])[1:nt]
        end
        ḡ = vec(sum(x .* hh, dims=2))          # Eq, 11 ĥ
        c = hh ./ ḡ .- 1                       # Eq. 12
        x .*= (1 .- c)                          # Eq. 16 (or 15?)
        # Convergence check
        err = maximum(abs, c)
        e3 = abs((e1 - err) / e1)
        e2 = e1
        e1 = err
    end
    return ḡ
end
function successive_flux(t, rb, xy, m::AbstractGroundModel; interp::Bool = true)
    @assert all(>(0), t) "t must be strictly positive to avoid a singular zero-time response"
    interp || return successive_flux(_response_array(t, rb, xy, m))
    return _subsampled_gfunction(t, tb -> successive_flux(_response_array(tb, rb, xy, m)))
end

"""
    bloc_matrix(g)
    bloc_matrix(t, rb, xy, m::AbstractGroundModel)

Function that computes the spatial superposition of a borefield using the bloc matrix approach
of Dusseault et al. (2018) to obtain g-functions of a borefield. This approach assumes that heat
flux is uniform along all the borehole, and that the mean temperature is the same for all
boreholes (Type II). The g-function generated is for an impulse of 1 W/m.
# Arguments
    - `g`: A 3D g-function matrix for all radius of the borefield (nt x nb x nb) [°Cm/W]
        - Each time step (1 x nb x nb) has the borefield response for each radius between boreholes.
        - Can be pre-computed with any ground model and passed directly.
    - `t`: Time vector (nt x 1) [s]
    - `rb`: Borehole radius [m]
    - `xy`: Matrix of borehole coordinates where the line source is at (0,0) (nb x 2) [m]
        - Can be computed with `borefield(:rectangle, ...)` from borefield.jl.
    - `m`: Ground model parameters (e.g. `FLSModel(150, 4, 3.0, 2e6)`)
# Output
    - `ḡ`: g-function of the borefield spatial superposition [°C·m/W]
# Reference
    - Dusseault, B., Pasquier, P., & Marcotte, D. (2018). A block matrix formulation for efficient
        g-function construction. Renewable Energy, 121, 249–260.
        https://doi.org/10.1016/j.renene.2017.12.092
"""
function bloc_matrix(gm::AbstractArray{<:Real})
    # Basic parameters
    nt, nb1, nb2 = size(gm)
    @assert nb1 == nb2 "g must be nt × nb × nb"
    @assert all(isfinite, gm) "g must contain only finite values"

    # Building the convolution matrix
    G = zeros(nt * nb1, nt * nb2)
    for i in 1:nt
        for ii in 1:nb1
            for jj in 1:nb2
                G[(ii-1)*nt+i:ii*nt, (jj-1)*nt+i] = gm[i:end, ii, jj]
            end
        end
    end

    # Create inputs to solve the linear system
    Gₕ = [[G; repeat(I(nt), 1, nb1)] [repeat(I(nt), nb2, 1); zeros(nt, nt)]]
    b = zeros((nb1 + 1) * nt)
    b[nb1*nt+1] = 1

    # Solve the linear system
    sol = Gₕ \ b
    # Output the transfer function
    return -sol[nb1*nt+1:end]
end
# Like `successive_flux`, the model-level method builds the response array through `_response_array`
# — the single per-model extension point — and feeds the low-level kernel. This is the instantaneous
# per-step formulation, so it needs no time sub-sampling of its own (`ground_response` still offers
# `interp` for it as a pure performance option).
bloc_matrix(t, rb, xy, m::AbstractGroundModel) = bloc_matrix(_response_array(t, rb, xy, m))

"""
    segment_response(g, Hseg, Href)
    segment_response(t, rb, xy, m::FLSModel)

Spatial superposition under the **segment boundary condition BC-III** of Cimmino and Bernier (2014):
each borehole is divided into `m.nseg` segments and the heat flux is allowed to vary from segment to
segment so that the **borehole-wall temperature is uniform over the whole field** (and constant in
depth). This is the finest of the three boundary conditions available in the package:

  - **BC-I** — equal heat flux on every borehole → [`uniform_flux`](@ref);
  - **BC-II** — uniform flux along each borehole, equal *mean* wall temperature →
    [`successive_flux`](@ref) / [`bloc_matrix`](@ref);
  - **BC-III** — flux varies *within* each borehole to keep a uniform wall temperature everywhere →
    `segment_response`.

The method generalises the block-matrix formulation of [`bloc_matrix`](@ref) from `nb` boreholes to
`nb·nseg` segments: it assembles one space–time block-Toeplitz system whose unknowns are the segment
fluxes and the common wall temperature, with an energy-balance row weighted by segment length. With
`nseg = 1` it reduces exactly to [`bloc_matrix`](@ref) (BC-II). As in the other methods, the internal
space–time solve produces the borefield g-function `ḡ(t)` for a unit step; load-side temporal
superposition is applied separately (by convolution) downstream.

# Arguments
    - `g`: A 4-argument segment g-array (nt × nS × nS) with `nS = nb·nseg`, `g[k, p, q]` the mean
        response over receiving segment `p` to a unit flux on emitting segment `q`. Segments are
        ordered borehole-major (`p = (borehole − 1)·nseg + segment`).
    - `Hseg`: Segment lengths (nS × 1) [m].
    - `Href`: Reference length used in the energy-balance constraint [m]; use the borehole length `H`
        so that `nseg = 1` reproduces [`bloc_matrix`](@ref).
    - `t`: Time vector (nt × 1) [s]
    - `rb`: Borehole radius [m]
    - `xy`: Matrix of borehole coordinates (nb × 2) [m] — can be built with `borefield(:rectangle, …)`.
    - `m`: Finite line source model with `m.nseg > 1` (e.g. `FLSModel(150, 4, 3.0, 2e6, 8)`).
# Output
    - `ḡ`: g-function of the borefield spatial superposition [°C·m/W]
# Reference
    - Cimmino, M., & Bernier, M. (2014). A semi-analytical method to generate g-functions for
        geothermal bore fields. International Journal of Heat and Mass Transfer, 70, 641–650.
        https://doi.org/10.1016/j.ijheatmasstransfer.2013.11.037
"""
function segment_response(g::AbstractArray{<:Real,3}, Hseg::AbstractVector{<:Real}, Href::Real)
    nt, nS1, nS2 = size(g)
    @assert nS1 == nS2 "g must be nt × nS × nS"
    @assert length(Hseg) == nS1 "Hseg must have one entry per segment"
    @assert all(isfinite, g) "g must contain only finite values"

    # Space–time block-Toeplitz convolution matrix (identical structure to bloc_matrix, but over
    # nS segments instead of nb boreholes).
    G = zeros(nt * nS1, nt * nS2)
    for i in 1:nt
        for ii in 1:nS1
            for jj in 1:nS2
                G[(ii-1)*nt+i:ii*nt, (jj-1)*nt+i] = g[i:end, ii, jj]
            end
        end
    end

    # Energy-balance row block: length-weighted sum of the segment (incremental) fluxes. Scaling the
    # weights by Hseg and the right-hand side by Href leaves the nseg = 1 case identical to
    # bloc_matrix (where every weight is 1 and the RHS is 1).
    W = repeat(Matrix{Float64}(I, nt, nt), 1, nS2)
    for jj in 1:nS2
        @views W[:, (jj-1)*nt+1:jj*nt] .*= Hseg[jj]
    end

    # Uniform-wall-temperature constraint (one common temperature for every segment).
    Gₕ = [[G; W] [repeat(I(nt), nS2, 1); zeros(nt, nt)]]
    b = zeros((nS1 + 1) * nt)
    b[nS1*nt+1] = Href

    sol = Gₕ \ b
    return -sol[nS1*nt+1:end]
end
function segment_response(t, rb, xy, m::FLSModel)
    g3D, Hseg, Href = _segment_g_array(t, rb, xy, m)
    return segment_response(g3D, Hseg, Href)
end

"""
    _segment_g_array(t, rb, xy, m::FLSModel)

Build the segment-to-segment g-array shared by both BC-III solvers. Each borehole is split into
`m.nseg` equal-length segments; the response of every (receiving, emitting) segment pair is assembled
into an `nt × nS × nS` array (`nS = nb·nseg`, borehole-major segment order). Returns the array along
with the segment lengths and the reference length used in the energy-balance constraint.
"""
function _segment_g_array(t, rb, xy, m::FLSModel)
    @assert all(>(0), t) "t must be strictly positive to avoid a singular zero-time response"
    nseg = m.nseg
    @assert nseg ≥ 1 "nseg must be ≥ 1"
    nb = size(xy, 1)
    H, D, ks, Cs = m.H, m.D, m.ks, m.Cs

    Hs   = H / nseg                                   # Segment length (equal-length segments)
    Dseg = [D + (k - 1) * Hs for k in 1:nseg]          # Buried depth of each segment top [m]
    nS   = nb * nseg

    # Horizontal borehole-to-borehole distances (diagonal = rb); segments of the same borehole share
    # the same horizontal position, so their separation is taken as rb (already on the diagonal).
    r, = borefield_geometry(xy, rb)
    rmat = max.(r, float(rb))

    # The segment response depends only on the triple (horizontal distance, emitting-segment depth,
    # receiving-segment depth). There are far fewer unique horizontal distances than borehole pairs,
    # so each (distance, source-depth, receiver-depth) response is evaluated once and scattered — the
    # same "similarities" reuse pygfunction relies on. For equal-length segments the depth pair is
    # symmetric (g[p,q] = g[q,p]), so only the js ≤ is triangle is computed.
    ru   = unique(round.(vec(rmat); digits = 10))
    ridx = Dict(ru[k] => k for k in eachindex(ru))
    nt   = length(t)
    cache = Array{Vector{Float64}}(undef, length(ru), nseg, nseg)
    for k in eachindex(ru)
        for is in 1:nseg, js in 1:is
            g = fls(t, ru[k], Hs, Dseg[js], Hs, Dseg[is], ks, Cs)
            cache[k, js, is] = g
            cache[k, is, js] = g
        end
    end

    # Segment index → (borehole, segment); borehole-major ordering.
    seg = [(bh, k) for bh in 1:nb for k in 1:nseg]

    g3D = zeros(nt, nS, nS)
    for q in 1:nS
        jb, js = seg[q]                                # emitting segment
        for p in 1:nS
            ib, is = seg[p]                            # receiving segment
            k = ridx[round(rmat[ib, jb]; digits = 10)]
            @views g3D[:, p, q] .= cache[k, js, is]
        end
    end

    return g3D, fill(Hs, nS), H
end

"""
    segment_response_marching(g, Hseg, Href)
    segment_response_marching(t, rb, xy, m::FLSModel; interp = true)

BC-III via **stepwise time-marching** (Cimmino & Bernier, 2014; Cimmino, 2018) — an alternative to
the block-matrix [`segment_response`](@ref). Instead of one large space–time system, it marches
forward one time step at a time, solving a small `(nS + 1) × (nS + 1)` system per step for the current
segment-flux increments and the common wall temperature. The coupling to the past is carried
explicitly by **temporal superposition** of the previously computed flux increments with the
segment-to-segment response (a discrete convolution, index-based on the given time grid, matching the
convention of [`successive_flux`](@ref)):

```math
T_i(t_p) = \\sum_{k=1}^{p} \\sum_j \\Delta q_j(t_k)\\, g_{ij}(t_{p-k+1}), \\qquad
\\sum_j H_j\\, q_j(t_p) = H_{ref}.
```

The two BC-III solvers use **different temporal formulations** — they are not merely two
discretisations of one problem:

  - [`segment_response`](@ref) (block matrix) applies the spatial superposition **independently at
    each time step** using ``g(t_p)`` as the kernel — it carries no load history (it is algebraically
    an instantaneous per-step solve). Cost `O((nS·nt)³)`, memory `O((nS·nt)²)`.
  - `segment_response_marching` applies **incremental temporal superposition** of the load history.
    Cost `O(nt·nS³ + nt²·nS²)`, memory `O(nS·nt)`.

They agree only when the flux never redistributes in time — i.e. `nseg = 1` on a geometrically
symmetric field. Otherwise (any `nseg > 1`, or an asymmetric field) they converge to **different**
limits a few percent apart; refining the time grid does not close the gap. The marching result is the
temporally rigorous one, and it is also far cheaper at scale. See
`dev/dev_bc3_marching_vs_block.jl` for a side-by-side comparison.

# Arguments
    - `g`: Segment g-array (nt × nS × nS), `g[k, p, q]` the response over receiving segment `p` to a
        unit flux on emitting segment `q`. Same layout as [`segment_response`](@ref).
    - `Hseg`: Segment lengths (nS × 1) [m].
    - `Href`: Reference length in the energy-balance constraint [m] (use the borehole length `H`).
    - `t`, `rb`, `xy`, `m`: as in [`segment_response`](@ref); `m::FLSModel` carries `nseg`.
    - `interp` (default `true`): like [`successive_flux`](@ref), this is a temporal solver whose
      index-based convolution needs a constant step, so `interp = true` re-grids onto the internal
      constant-step sub-sampling and interpolates to `t` — a correctness requirement on non-uniform
      `t`. Pass `interp = false` to solve directly on a uniform `t`.
# Output
    - `ḡ`: g-function of the borefield spatial superposition [°C·m/W]
# Reference
    - Cimmino, M. (2018). Fast calculation of the g-functions of geothermal borehole fields using
        similarities in the evaluation of the finite line source solution. Journal of Building
        Performance Simulation, 11(6), 655–668. https://doi.org/10.1080/19401493.2017.1423390
    - Cimmino, M., & Bernier, M. (2014). A semi-analytical method to generate g-functions for
        geothermal bore fields. International Journal of Heat and Mass Transfer, 70, 641–650.
"""
function segment_response_marching(g::AbstractArray{<:Real,3}, Hseg::AbstractVector{<:Real},
    Href::Real)
    nt, nS1, nS2 = size(g)
    @assert nS1 == nS2 "g must be nt × nS × nS"
    @assert length(Hseg) == nS1 "Hseg must have one entry per segment"
    @assert all(isfinite, g) "g must contain only finite values"
    T = float(eltype(g))
    nS = nS1

    # Per-step system  [ g(t₁)  −1 ; Hseg  0 ] · [Δq; T] = [−T₀; (Href at p=1 else 0)].
    # The matrix is constant across steps (the current-step response is g at the first grid lag), so
    # it is factorised once and reused; only the right-hand side (the load-history term) changes.
    A = zeros(T, nS + 1, nS + 1)
    @views A[1:nS, 1:nS] .= g[1, :, :]
    @views A[1:nS, nS+1] .= -one(T)
    @views A[nS+1, 1:nS] .= Hseg
    F = lu(A)

    dQ  = zeros(T, nt, nS)        # segment-flux increments per step
    gf  = zeros(T, nt)           # field g-function = common wall temperature
    T0  = zeros(T, nS)
    rhs = zeros(T, nS + 1)
    for p in 1:nt
        # Temporal superposition of the past flux increments with the response at each elapsed lag.
        fill!(T0, zero(T))
        for k in 1:p-1
            gl = @view g[p-k+1, :, :]
            dq = @view dQ[k, :]
            mul!(T0, gl, dq, one(T), one(T))   # T0 .+= gl * dq
        end
        @views rhs[1:nS] .= .-T0
        rhs[nS+1] = p == 1 ? T(Href) : zero(T)
        x = F \ rhs
        @views dQ[p, :] .= x[1:nS]
        gf[p] = x[nS+1]
    end
    return gf
end
# Like `successive_flux`, the marching solver is temporal (index-based convolution) and therefore
# assumes a constant time step; `interp = true` (default) solves on the internal constant-step grid
# and interpolates to `t` — a correctness requirement on any non-uniform `t`, not just a speed-up.
# Pass `interp = false` to solve directly on `t` (uniform grids only).
function segment_response_marching(t, rb, xy, m::FLSModel; interp::Bool = true)
    @assert all(>(0), t) "t must be strictly positive to avoid a singular zero-time response"
    solve(tb) = segment_response_marching(_segment_g_array(tb, rb, xy, m)...)
    interp || return solve(t)
    return _subsampled_gfunction(t, solve)
end