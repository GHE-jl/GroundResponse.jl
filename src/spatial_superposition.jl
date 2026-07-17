using LinearAlgebra
using DSP: conv

"""
    uniform_flux(g)
    uniform_flux(t, rb, xy, m::AbstractGroundModel)

Spatial superposition under a **uniform, equal heat flux on every borehole** (the M1 model of
Guo et al., 2021), equivalent to boundary condition Type-I (BC-I).
Every borehole emits the same constant flux; the response of borehole `i` is the sum of its
self response and the mutual responses from all other boreholes, and the borefield g-function is
the average of these responses over all boreholes (Guo et al., 2021, Eqs. 11–12).
Contrary to `successive_flux` / `bloc_matrix`, which enforce a uniform mean borehole-wall
temperature (BC-II) through a linear solve, M1 requires no matrix inversion. It is therefore
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
function uniform_flux(t, rb::Real, xy::AbstractMatrix{<:Real}, m::AbstractGroundModel)
    _check_time(t)
    return uniform_flux(_borehole_response(t, rb, xy, m))
end

"""
    successive_flux(g)
    successive_flux(t, rb, xy, m::AbstractGroundModel; interp = true)

Iteratively solve spatial superposition for a borefield using the successive flux approach of
Nguyen and Pasquier (2021) to obtain the g-functions of a borefield. This approach assumes that heat
flux is uniform along all the borehole, and that the mean temperature is the same for all
boreholes (BC-II). The g-function generated is for an impulse of 1 W/m.

This is a **temporal** solver: it reproduces the load history by a spectral convolution, which
assumes a **constant time step**. The `interp` keyword (default `true`) solves on an internal
constant-step grid and PCHIP-interpolates to the requested `t`. For this solver `interp = true` is a
**correctness** requirement on any non-uniform (e.g. log-spaced) `t`. Pass `interp = false` only 
when `t` is already uniformly spaced, to solve directly on it.
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
function successive_flux(t, rb, xy::AbstractMatrix{<:Real}, m::AbstractGroundModel;
    interp::Bool = true)
    _check_time(t)
    interp || return successive_flux(_borehole_response(t, rb, xy, m))
    return _interp_gfunction(t, tb -> successive_flux(_borehole_response(tb, rb, xy, m)))
end

"""
    bloc_matrix(g)
    bloc_matrix(t, rb, xy, m::AbstractGroundModel)

Function that computes the spatial superposition of a borefield using the bloc matrix approach
of Dusseault et al. (2018) to obtain g-functions of a borefield. This approach assumes that heat
flux is uniform along all the borehole, and that the mean temperature is the same for all
boreholes (BC-II). The g-function generated is for an impulse of 1 W/m.
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
function bloc_matrix(t, rb, xy::AbstractMatrix{<:Real}, m::AbstractGroundModel)
    _check_time(t)
    return bloc_matrix(_borehole_response(t, rb, xy, m))
end

"""
    _segment_g(t, rb, xy, m::FLSModel)

Build the segment-to-segment g-array shared by both BC-III solvers, the FLS-only counterpart of
[`borehole_response`](@ref). Each borehole is split into `m.nseg` equal-length segments; th
response of every (receiving, emitting) segment pair is assembled into an `nt × nS × nS` array
(`nS = nb·nseg`, borehole-major segment order). Returns the array along with the segment lengths and
the reference length used in the energy-balance constraint.
"""
function _segment_g(t, rb, xy::AbstractMatrix{<:Real}, m::FLSModel)
    nseg = m.nseg
    @assert nseg ≥ 1 "nseg must be ≥ 1"
    nb = size(xy, 1)
    H, D, ks, Cs = m.H, m.D, m.ks, m.Cs

    Hs   = H / nseg                                   # Segment length (equal-length segments)
    Dseg = [D + (k - 1) * Hs for k in 1:nseg]          # Buried depth of each segment top [m]
    nS   = nb * nseg

    # Horizontal borehole-to-borehole distances
    r, = borefield_geometry(xy, rb)
    rmat = max.(r, float(rb))

    # Precompute the segment-to-segment g-functions for every unique borehole distance
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

    # Segment index → (borehole, segment)
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
    segment_response(g, Hseg, Href)
    segment_response(t, rb, xy, m::FLSModel)

Spatial superposition under the **segment boundary condition BC-III** of Cimmino and Bernier (2014):
each borehole is divided into `m.nseg` segments and the heat flux is allowed to vary from segment to
segment so that the **borehole-wall temperature is uniform over the whole field** (and constant in
depth).
The method generalises the block-matrix formulation of [`bloc_matrix`](@ref) from `nb` boreholes to
`nb·nseg` segments: it assembles one space–time block-Toeplitz system whose unknowns are the segment
fluxes and the common wall temperature, with an energy-balance row weighted by segment length. With
`nseg = 1` it reduces exactly to [`bloc_matrix`](@ref) (BC-II).
# Arguments
    - `g`: A 4-argument segment g-array (nt × nS × nS) with `nS = nb·nseg`, `g[k, p, q]` the mean
        response over receiving segment `p` to a unit flux on emitting segment `q`. Segments are
        ordered borehole-major (`p = (borehole − 1)·nseg + segment`).
    - `Hseg`: Segment lengths (nS × 1) [m].
    - `Href`: Reference length used in the energy-balance constraint [m]; use the borehole length
        `H` so that `nseg = 1` reproduces [`bloc_matrix`](@ref).
    - `t`: Time vector (nt × 1) [s]
    - `rb`: Borehole radius [m]
    - `xy`: Matrix of borehole coordinates where the line source is at (0,0) (nb x 2) [m]
        - Can be computed with `borefield(:rectangle, ...)` from borefield.jl.
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

    # Space–time block-Toeplitz convolution matrix
    G = zeros(nt * nS1, nt * nS2)
    for i in 1:nt
        for ii in 1:nS1
            for jj in 1:nS2
                G[(ii-1)*nt+i:ii*nt, (jj-1)*nt+i] = g[i:end, ii, jj]
            end
        end
    end

    # Energy-balance row block: length-weighted sum of the segment (incremental) fluxes.
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
function segment_response(t, rb, xy::AbstractMatrix{<:Real}, m::FLSModel)
    _check_time(t)
    g3D, Hseg, Href = _segment_g(t, rb, xy, m)
    return segment_response(g3D, Hseg, Href)
end

"""
    segment_response_marching(g, Hseg, Href)
    segment_response_marching(t, rb, xy, m::FLSModel; interp = true)

BC-III via **stepwise time-marching** (Cimmino & Bernier, 2014; Cimmino, 2018), an alternative to
the block-matrix [`segment_response`](@ref). Instead of one large space–time system, it marches
forward one time step at a time, solving a small `(nS + 1) × (nS + 1)` system per step for the
current segment-flux increments and the common wall temperature.
The two BC-III solvers use **different temporal formulations**:
  - [`segment_response`](@ref) (block matrix) applies the spatial superposition **independently at
    each time step** using ``g(t_p)`` as the kernel.
  - `segment_response_marching` applies **incremental temporal superposition** of the load history.
They agree only when the flux never redistributes in time — i.e. `nseg = 1` on a geometrically
symmetric field. Otherwise (any `nseg > 1`, or an asymmetric field) they converge to **different**
limits a few percent apart; refining the time grid does not close the gap.
# Arguments
    - `g`: Segment g-array (nt × nS × nS), `g[k, p, q]` the response over receiving segment `p` to a
        unit flux on emitting segment `q`.
    - `Hseg`: Segment lengths (nS × 1) [m].
    - `Href`: Reference length in the energy-balance constraint [m] (use the borehole length `H`).
    - `t`: Time vector (nt × 1) [s]
    - `rb`: Borehole radius [m]
    - `xy`: Matrix of borehole coordinates where the line source is at (0,0) (nb x 2) [m]
        - Can be computed with `borefield(:rectangle, ...)` from borefield.jl.
    - `m`: Finite line source model with `m.nseg > 1` (e.g. `FLSModel(150, 4, 3.0, 2e6, 8)`).
    - `interp` (default `true`): A correctness requirement on non-uniform `t`. Pass `interp = false`
        to solve directly on a uniform `t`.
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

    # Per-step system
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
function segment_response_marching(t, rb, xy::AbstractMatrix{<:Real}, m::FLSModel; interp::Bool = true)
    _check_time(t)
    solve(tb) = segment_response_marching(_segment_g(tb, rb, xy, m)...)
    interp || return solve(t)
    return _interp_gfunction(t, solve)
end