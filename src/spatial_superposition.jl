using LinearAlgebra
using DSP: conv

# Spatial superposition of the borehole responses into a borefield g-function, by boundary
# condition:
#   - `bc_I`: equal, uniform heat flux on every borehole (no solve);
#   - `bc_II`: uniform flux along each borehole, equal mean wall temperature;
#   - `bc_III`: flux free to vary between the segments of each borehole, uniform wall temperature.
# BC-II and BC-III share one direct solver, `_marching_solver`. The previous solvers are kept at the
# end of the file under a `_legacy` suffix.

"""
    _marching_solver(g, Hseg, Href)

Time-marching solution of the spatio-temporal superposition with a common wall temperature, shared
by [`bc_II`](@ref) (one source per borehole) and [`bc_III`](@ref) (one source per segment).

The flux of every source is constant over each step. The unknowns of step `p` are the flux
increments `f_p = q(t_p) − q(t_{p−1})` and the common wall temperature `ḡ_p`. With
`g(t_p − t_{l−1}) = g[p−l+1]` on a uniform grid, the temperature of every source is the sum of the
current step and of the load history:

    ḡ_p = g[1] ⋅ f_p + h_p,        h_p = Σ_{l<p} g[p−l+1] ⋅ f_l          (known)

and each step solves the small `(nS + 1) × (nS + 1)` system

    [ g[1]  −1 ] [ f_p ]   [     −h_p      ]
    [ Hsegᵀ  0 ] [ ḡ_p ] = [ Href if p = 1 ]   (0 otherwise: the total heat rate stays constant)

This is the block-matrix system of Dusseault et al. (2018) solved by block forward substitution, so
the solution is the same, without iteration, at a fraction of the cost and memory. `g[1]` is the
same at every step, so it is factorised once.

The forward substitution amplifies rounding errors by `g[k]/g[1]` at every step, so it is unstable
when the first step lies in the onset of the self response (`Δt ≲ rb²/α`, below about one hour for
usual boreholes), where `g` grows by orders of magnitude between the first steps. Heat has not
left the borehole wall, the sources do not interact and their fluxes stay uniform. Then
when a self response grows more than tenfold over the first step, the uniform-flux response (as
[`bc_I`](@ref)) is returned instead.
# Arguments
    - `g`: Response array (nt × nS × nS) on a uniform time grid `t = Δt, 2Δt, …, ntΔt`, `g[k, p, q]`
        the response of source `p` to a unit flux on source `q`. May be asymmetric.
    - `Hseg`: Length of each source (nS × 1) [m]
    - `Href`: Total length the heat rate is referred to [m]
# Output
    - `ḡ`: Common wall g-function at each time (nt × 1) for a total heat rate `Href` [°C·m/W]
"""
function _marching_solver(g::AbstractArray{<:Real,3}, Hseg::AbstractVector{<:Real}, Href::Real)
    # Initiation
    nt, nS1, nS2 = size(g)
    @assert nS1 == nS2 "g must be nt × nS × nS"
    nS = nS1
    @assert length(Hseg) == nS "Hseg must have one entry per source"
    @assert all(isfinite, g) "g must contain only finite values"
    T = float(eltype(g))

    # Contiguous lag slices: G[:, :, k] = g[k, :, :]
    G = permutedims(g, (2, 3, 1))

    # Onset of the self response: no interaction, uniform flux (the marching would be unstable)
    if nt > 1 && any(i -> G[i, i, 2] > 10 * G[i, i, 1], 1:nS)
        q = Href / sum(Hseg)
        return [q * dot(Hseg, sum(view(G, :, :, k); dims = 2)) / sum(Hseg) for k in 1:nt]
    end

    # Step system, the same at every step: factorised once
    A = zeros(T, nS + 1, nS + 1)
    @views A[1:nS, 1:nS] .= G[:, :, 1]
    @views A[1:nS, nS+1] .= -one(T)
    @views A[nS+1, 1:nS] .= Hseg
    F = lu!(A)

    # Marching loop: solve for flux increments and common wall g-function at each time step
    f   = zeros(T, nS, nt)          # flux increments
    h   = zeros(T, nS)              # load history
    x   = zeros(T, nS + 1)          # right-hand side, then solution of the step
    out = zeros(T, nt)
    for p in 1:nt
        fill!(h, zero(T))
        for l in 1:p-1
            mul!(h, view(G, :, :, p - l + 1), view(f, :, l), one(T), one(T))
        end
        @views x[1:nS] .= .-h
        x[nS+1] = p == 1 ? T(Href) : zero(T)
        ldiv!(F, x)
        @views f[:, p] .= x[1:nS]
        out[p] = x[nS+1]
    end
    return out
end

"""
    bc_I(g)
    bc_I(t, rb, xy, m::AbstractGroundModel; interp = true)

Spatial superposition under boundary condition **BC-I**: every borehole emits the same constant flux
(the M1 model of Guo et al., 2021). The response of borehole `i` is the sum of its self response and
of the responses to all other boreholes, and the borefield g-function is their average (Guo et al.,
2021, Eqs. 11–12). No system is solved, so it is robust for the asymmetric responses of the
advection models ([`MILSModel`](@ref), [`MFLSModel`](@ref)).
# Arguments
    - `g`: Response array (nt × nb × nb), or matrix (nb × nb) for a single time. `g[k, i, j]` is the
        response at borehole `i` to a unit flux at borehole `j`. May be asymmetric.
    - `t`: Time vector (nt × 1) [s]
    - `rb`: Borehole radius [m]
    - `xy`: Borehole coordinates (nb × 2) [m], e.g. from [`borefield`](@ref)
    - `m`: Ground model (e.g. `FLSModel(150, 4, 3.0, 2e6)`)
# Keywords
    - `interp` (default `true`): compute on the internal sub-sampling grid and PCHIP-interpolate to
        `t`. A performance choice only for BC-I, which is exact at any `t`.
# Output
    - `g`: Borefield g-function for 1 W/m on every borehole [°C·m/W]
# Reference
    - Guo, Y., Hu, X., Banks, J., & Liu, W. V. (2021). Considering buried depth for vertical
        borehole heat exchangers in a borehole field with groundwater flow — An extended solution.
        Energy and Buildings, 235, 110722. https://doi.org/10.1016/j.enbuild.2021.110722
"""
function bc_I(g::AbstractMatrix{<:Real})
    # 1 time step, g is nb × nb
    nb1, nb2 = size(g)
    @assert nb1 == nb2 "g must be nb × nb"
    @assert all(isfinite, g) "g must contain only finite values"
    return sum(g) / nb1
end
function bc_I(g::AbstractArray{<:Real,3})
    # Multiple time steps, g is nt × nb × nb
    nt, nb1, nb2 = size(g)
    @assert nb1 == nb2 "g must be nt × nb × nb"
    @assert all(isfinite, g) "g must contain only finite values"
    return [sum(@view g[k, :, :]) / nb1 for k in 1:nt]
end
function bc_I(t, rb, xy::AbstractMatrix{<:Real}, m::AbstractGroundModel; interp::Bool = true)
    # Single or multiple time steps: compute the borehole responses and apply BC-I superposition
    _check_time(t)
    return _apply_interp(t, interp, tt -> bc_I(_borehole_response(tt, rb, xy, m))) #TODO this is the only place where _apply_interp is used. Maybe simplify and only use the other function?
end

"""
    bc_II(g)
    bc_II(t, rb, xy, m::AbstractGroundModel; interp = true)

Spatial superposition under boundary condition **BC-II**: the flux is uniform along each borehole,
the mean wall temperature is the same for every borehole, and the total heat rate is 1 W/m per
borehole. The flux of each borehole is free to change in time, so the load history is superposed in
time and in space. Solved directly by time marching (`_marching_solver`), with the same
solution as the block matrix of Dusseault et al. (2018).
# Arguments
    - `g`: Response array (nt × nb × nb) on a uniform time grid `t = Δt, 2Δt, …, ntΔt`. `g[k, i, j]`
        is the response at borehole `i` to a unit flux at borehole `j`. May be asymmetric.
    - `t`: Time vector (nt × 1) [s]
    - `rb`: Borehole radius [m]
    - `xy`: Borehole coordinates (nb × 2) [m], e.g. from [`borefield`](@ref)
    - `m`: Ground model (e.g. `FLSModel(150, 4, 3.0, 2e6)`)
# Keywords
    - `interp` (default `true`): solve on the internal sub-sampling grid (uniform blocks, see
        `time_sampling.jl`) and PCHIP-interpolate to `t`. Required for any `t` that is not
        `Δt, 2Δt, …, nΔt`; `interp = false` solves directly on such a `t`.
# Output
    - `g`: Borefield g-function for 1 W/m on every borehole [°C·m/W]
# References
    - Dusseault, B., Pasquier, P., & Marcotte, D. (2018). A block matrix formulation for efficient
        g-function construction. Renewable Energy, 121, 249–260.
        https://doi.org/10.1016/j.renene.2017.12.092
    - Cimmino, M., & Bernier, M. (2014). A semi-analytical method to generate g-functions for
        geothermal bore fields. International Journal of Heat and Mass Transfer, 70, 641–650.
        https://doi.org/10.1016/j.ijheatmasstransfer.2013.11.037
"""
function bc_II(g::AbstractArray{<:Real,3})
    # Single time step: g is nb × nb
    nb = size(g, 2)
    return _marching_solver(g, ones(nb), nb)
end
function bc_II(t, rb, xy::AbstractMatrix{<:Real}, m::AbstractGroundModel; interp::Bool = true)
    # Multiple time steps: g is nt × nb × nb
    _check_time(t)
    solve(tt) = bc_II(_borehole_response(tt, rb, xy, m))
    interp || return (_check_uniform_steps(t); solve(t))
    return _interp_gfunction(t, solve)
end

"""
    bc_III(g, Hseg, Href)
    bc_III(t, rb, xy, m::FLSModel; interp = true)

Spatial superposition under boundary condition **BC-III** (Cimmino & Bernier, 2014): each borehole
is split into `m.nseg` equal segments whose fluxes are free to vary, so that the wall temperature is
uniform over the whole field. The same marching solver as [`bc_II`](@ref), over `nb·nseg` segments,
with an energy balance weighted by the segment lengths. With `nseg = 1` it is identical to
[`bc_II`](@ref). This BC can only be applied to finite length ground models.
# Arguments
    - `g`: Segment response array (nt × nS × nS), `nS = nb·nseg`, on a uniform time grid
        `t = Δt, 2Δt, …, ntΔt`. `g[k, p, q]` is the mean response over segment `p` to a unit flux on
        segment `q`, segments ordered borehole-major (`p = (borehole − 1)·nseg + segment`).
        Obtained from `_segment_g`.
    - `Hseg`: Segment lengths (nS × 1) [m]
    - `Href`: Reference length of the energy balance [m]: the total length `nb·H` gives 1 W/m per
        borehole.
    - `t`: Time vector (nt × 1) [s]
    - `rb`: Borehole radius [m]
    - `xy`: Borehole coordinates (nb × 2) [m], e.g. from [`borefield`](@ref)
    - `m`: Finite line source model with its number of segments (e.g. `FLSModel(150, 4, 3.0, 2e6, 8)`)
# Keywords
    - `interp` (default `true`): as for [`bc_II`](@ref).
# Output
    - `g`: Borefield g-function for 1 W/m on every borehole [°C·m/W]
# Reference
    - Cimmino, M., & Bernier, M. (2014). A semi-analytical method to generate g-functions for
        geothermal bore fields. International Journal of Heat and Mass Transfer, 70, 641–650.
        https://doi.org/10.1016/j.ijheatmasstransfer.2013.11.037
"""
function bc_III(g::AbstractArray{<:Real,3}, Hseg::AbstractVector{<:Real}, Href::Real)
    # Marching solver for BC-III with segment response array, segment lengths, and reference length
    _marching_solver(g, Hseg, Href)
end
function bc_III(t, rb, xy::AbstractMatrix{<:Real}, m::FLSModel; interp::Bool = true)
    # Marching solver for BC-III with FLS model, borehole coordinates, and time vector
    _check_time(t)
    solve(tt) = bc_III(_segment_g(tt, rb, xy, m)...)
    interp || return (_check_uniform_steps(t); solve(t))
    return _interp_gfunction(t, solve)
end

"""
    _segment_g(t, rb, xy, m::FLSModel)

Build the segment-to-segment response array of BC-III, the FLS-only counterpart of
`_borehole_response`. Each borehole is split into `m.nseg` equal-length segments; the response of
every (receiving, emitting) segment pair is assembled into an `nt × nS × nS` array (`nS = nb·nseg`,
borehole-major segment order). Returns the array along with the segment lengths and the reference
length used in the energy-balance constraint.
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

    return g3D, fill(Hs, nS), nb * H
end

# Legacy solvers, superseded by `bc_II` / `bc_III` (same boundary conditions). Kept for comparison;
# the names without the `_legacy` suffix are deprecated aliases (see deprecated.jl).

"""
    successive_flux_legacy(g)
    successive_flux_legacy(t, rb, xy, m::AbstractGroundModel; interp = true)

Legacy BC-II solver, superseded by [`bc_II`](@ref). Iterative successive flux estimation of Nguyen
and Pasquier (2021): the fluxes of every borehole at every time are corrected until the wall
temperatures agree within 0.1 % (`max |T_i/T̄ − 1| < 1e-3`). The initial guess solves each time step
without the load history; the load history is then superposed by one convolution per borehole pair
(`DSP.conv`). It converges towards the [`bc_II`](@ref) solution, within about 1e-5 to 3e-4.

The convolution assumes a **constant time step** `t = Δt, 2Δt, …, nΔt`; `interp = true` (default)
provides it for any `t`. Pass `interp = false` only for such a `t`.
# Arguments
    - `g`: Response array (nt × nb × nb) [°C·m/W]
    - `t`: Time vector (nt × 1) [s]
    - `rb`: Borehole radius [m]
    - `xy`: Borehole coordinates (nb × 2) [m]
    - `m`: Ground model (e.g. `FLSModel(150, 4, 3.0, 2e6)`)
# Output
    - `g`: Borefield g-function for 1 W/m on every borehole [°C·m/W]
# Reference
    - Nguyen, A., & Pasquier, P. (2021). A successive flux estimation method for rapid g-function
        construction of small to large-scale ground heat exchanger. Renewable Energy, 165, 359–368.
        https://doi.org/10.1016/j.renene.2020.10.074
"""
function successive_flux_legacy(g::AbstractArray{<:Real,3})
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
    ḡ = zeros(eltype(g), nt)

    while e3 > 0.15 && e1 > 1e-3 && e1 < e2 && k < kmax
        k += 1                                  # Iteration counter
        f = vcat(x[1:1, :], diff(x, dims=1))    # Step fluxes (Eq. 4)
        # Temperature response via pairwise convolutions (Eq. 6, 8)
        hh = zeros(eltype(g), nt, nb)
        for j in 1:nb, i in 1:nb
            hh[:, i] .+= conv(f[:, j], g[:, i, j])[1:nt]
        end
        ḡ = vec(sum(x .* hh, dims=2))          # Eq, 11 ĥ
        c = hh ./ ḡ .- 1                       # Eq. 12
        x .*= (1 .- c)                          # Eq. 16 (or 15?)
        # Convergence check
        err = maximum(abs, c)
        e3 = abs((e1 - err) / e1)
        e2 = e1
        e1 = err
    end
    return nb .* ḡ
end
function successive_flux_legacy(t, rb, xy::AbstractMatrix{<:Real}, m::AbstractGroundModel;
    interp::Bool = true)
    _check_time(t)
    interp || return successive_flux_legacy(_borehole_response(t, rb, xy, m))
    return _interp_gfunction(t, tb -> successive_flux_legacy(_borehole_response(tb, rb, xy, m)))
end

"""
    bloc_matrix_legacy(g)
    bloc_matrix_legacy(t, rb, xy, m::AbstractGroundModel; interp = true)

Legacy BC-II solver, superseded by [`bc_II`](@ref) (same solution). Block-matrix formulation of
Dusseault et al. (2018): one linear system of size `(nb + 1)·nt` holds the temporal and spatial
superposition of every borehole at every time step, with a lower-triangular Toeplitz block
`G[k, l] = g(t_k − t_{l−1}) = g[k − l + 1]` per pair of boreholes. Its cost grows as `(nb·nt)³`.

The Toeplitz block requires a **constant time step** `t = Δt, 2Δt, …, nΔt`; `interp = true`
(default) provides it for any `t`. Pass `interp = false` only for such a `t`.
# Arguments
    - `g`: Response array (nt × nb × nb) [°C·m/W]
    - `t`: Time vector (nt × 1) [s]
    - `rb`: Borehole radius [m]
    - `xy`: Borehole coordinates (nb × 2) [m]
    - `m`: Ground model (e.g. `FLSModel(150, 4, 3.0, 2e6)`)
# Output
    - `g`: Borefield g-function for 1 W/m on every borehole [°C·m/W]
# Reference
    - Dusseault, B., Pasquier, P., & Marcotte, D. (2018). A block matrix formulation for efficient
        g-function construction. Renewable Energy, 121, 249–260.
        https://doi.org/10.1016/j.renene.2017.12.092
"""
function bloc_matrix_legacy(gm::AbstractArray{<:Real})
    # Basic parameters
    nt, nb1, nb2 = size(gm)
    @assert nb1 == nb2 "g must be nt × nb × nb"
    @assert all(isfinite, gm) "g must contain only finite values"

    # Building the convolution matrix: the flux increment of step i acts for t_k − t_{i−1}
    G = zeros(nt * nb1, nt * nb2)
    for i in 1:nt
        for ii in 1:nb1
            for jj in 1:nb2
                G[(ii-1)*nt+i:ii*nt, (jj-1)*nt+i] = gm[1:nt-i+1, ii, jj]
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
    return -nb1 .* sol[nb1*nt+1:end]
end
function bloc_matrix_legacy(t, rb, xy::AbstractMatrix{<:Real}, m::AbstractGroundModel;
    interp::Bool = true)
    _check_time(t)
    interp || return bloc_matrix_legacy(_borehole_response(t, rb, xy, m))
    return _interp_gfunction(t, tb -> bloc_matrix_legacy(_borehole_response(tb, rb, xy, m)))
end

"""
    segment_response_legacy(g, Hseg, Href)
    segment_response_legacy(t, rb, xy, m::FLSModel; interp = true)

Legacy BC-III solver, superseded by [`bc_III`](@ref) (same solution). The block matrix of
[`bloc_matrix_legacy`](@ref) generalised from `nb` boreholes to `nb·nseg` segments (Cimmino &
Bernier, 2014), with an energy balance weighted by the segment lengths. Its cost grows as
`(nb·nseg·nt)³`.

The Toeplitz time block requires a **constant time step** `t = Δt, 2Δt, …, nΔt`; `interp = true`
(default) provides it for any `t`. Pass `interp = false` only for such a `t`.
# Arguments
    - `g`: Segment response array (nt × nS × nS), `nS = nb·nseg`, segments ordered borehole-major
    - `Hseg`: Segment lengths (nS × 1) [m]
    - `Href`: Reference length of the energy balance [m] (the total length `nb·H`)
    - `t`: Time vector (nt × 1) [s]
    - `rb`: Borehole radius [m]
    - `xy`: Borehole coordinates (nb × 2) [m]
    - `m`: Finite line source model with its number of segments (e.g. `FLSModel(150, 4, 3.0, 2e6, 8)`)
# Output
    - `g`: Borefield g-function for 1 W/m on every borehole [°C·m/W]
# Reference
    - Cimmino, M., & Bernier, M. (2014). A semi-analytical method to generate g-functions for
        geothermal bore fields. International Journal of Heat and Mass Transfer, 70, 641–650.
        https://doi.org/10.1016/j.ijheatmasstransfer.2013.11.037
"""
function segment_response_legacy(g::AbstractArray{<:Real,3}, Hseg::AbstractVector{<:Real},
    Href::Real)
    nt, nS1, nS2 = size(g)
    @assert nS1 == nS2 "g must be nt × nS × nS"
    @assert length(Hseg) == nS1 "Hseg must have one entry per segment"
    @assert all(isfinite, g) "g must contain only finite values"

    # Space–time block-Toeplitz convolution matrix
    G = zeros(nt * nS1, nt * nS2)
    for i in 1:nt
        for ii in 1:nS1
            for jj in 1:nS2
                G[(ii-1)*nt+i:ii*nt, (jj-1)*nt+i] = g[1:nt-i+1, ii, jj]
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
function segment_response_legacy(t, rb, xy::AbstractMatrix{<:Real}, m::FLSModel;
    interp::Bool = true)
    _check_time(t)
    solve(tb) = segment_response_legacy(_segment_g(tb, rb, xy, m)...)
    interp || return solve(t)
    return _interp_gfunction(t, solve)
end

"""
    segment_response_marching_legacy(g, Hseg, Href)
    segment_response_marching_legacy(t, rb, xy, m::FLSModel; interp = true)

Legacy BC-III time-marching solver, superseded by [`bc_III`](@ref): the same algorithm, but each
lag of the load history is read as a strided slice `g[k, :, :]`, so it is slower on large problems.

It requires a **constant time step** `t = Δt, 2Δt, …, nΔt`; `interp = true` (default) provides it
for any `t`. Pass `interp = false` only for such a `t`.
# Arguments
    - `g`: Segment response array (nt × nS × nS), `nS = nb·nseg`, segments ordered borehole-major
    - `Hseg`: Segment lengths (nS × 1) [m]
    - `Href`: Reference length of the energy balance [m] (the total length `nb·H`)
    - `t`: Time vector (nt × 1) [s]
    - `rb`: Borehole radius [m]
    - `xy`: Borehole coordinates (nb × 2) [m]
    - `m`: Finite line source model with its number of segments (e.g. `FLSModel(150, 4, 3.0, 2e6, 8)`)
# Output
    - `g`: Borefield g-function for 1 W/m on every borehole [°C·m/W]
# Reference
    - Cimmino, M. (2018). Fast calculation of the g-functions of geothermal borehole fields using
        similarities in the evaluation of the finite line source solution. Journal of Building
        Performance Simulation, 11(6), 655–668. https://doi.org/10.1080/19401493.2017.1423390
"""
function segment_response_marching_legacy(g::AbstractArray{<:Real,3}, Hseg::AbstractVector{<:Real},
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
function segment_response_marching_legacy(t, rb, xy::AbstractMatrix{<:Real}, m::FLSModel;
    interp::Bool = true)
    _check_time(t)
    solve(tb) = segment_response_marching_legacy(_segment_g(tb, rb, xy, m)...)
    interp || return solve(t)
    return _interp_gfunction(t, solve)
end
