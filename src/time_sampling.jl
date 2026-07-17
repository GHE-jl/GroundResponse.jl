# The constant-step sub-sampling grid and PCHIP interpolation used to (a) satisfy the constant-step
# requirement of the temporal spatial-superposition solvers and (b) bound the cost of any backend on
# a large `t`.

# Fixed sub-sampling scheme for the TEMPORAL spatial-superposition solvers (`successive_flux`,
# `segment_response_marching`). Their spectral / index-based convolution assumes a constant time
# step, so passing an arbitrary (e.g. log-spaced) `t` silently corrupts the temporal superposition.
# The g-function is instead solved on the blocks below, each UNIFORMLY spaced, so the convolution
# within a block is valid, following the geometric hour→decade progression of Nguyen & Pasquier
# (2021).
# The same routine also serves the instantaneous / direct paths (`uniform_flux`, `bloc_matrix`,
# single-borehole models), where it is a pure performance approximation rather than a correctness
# requirement.

const _SUBSAMPLE_BLOCKS = ((2.0,          30),   # every 2 s   → 60 s
                           (120.0,        30),   # every 2 min → 1 h
                           (3600.0,       24),   # every 1 h   → 1 day
                           (86_400.0,     30),   # every 1 day → 1 month (30-day)
                           (2_592_000.0,  12),   # every 1 month → 1 year (360-day)
                           (31_104_000.0, 10))   # every 1 year → 1 decade

const _SUBSAMPLE_DECADE = 311_040_000.0          # 1 decade [s]; count grows to cover the request

"""
    _append_strictly_increasing!(ts, gs, tb, gb)

Returns `ts` and `gs` with the contents of `tb` and `gb` appended, skipping any entries in `tb`
that are not strictly greater than the last entry of `ts`. This is used to build a strictly
increasing time vector for the sub-sampled g-function solver, which is required for the PCHIP
interpolation to work correctly.
# Arguments
    - `ts`: Vector of time steps [s]
    - `gs`: Vector of g-function values [°Cm/W]
    - `tb`: Vector of new time steps to append [s]
    - `gb`: Vector of new g-function values to append [°Cm/W]
# Output
    - `ts`: Updated vector of time steps [s]
    - `gs`: Updated vector of g-function values [°Cm/W]
"""
function _append_strictly_increasing!(ts, gs, tb, gb)
    @inbounds for k in eachindex(tb)
        if isempty(ts) || tb[k] > ts[end] * (1 + 1e-12)
            push!(ts, tb[k])
            push!(gs, gb[k])
        end
    end
    return nothing
end

"""
    _interp_gfunction(t, solve)

Solve a temporal g-function on the fixed sub-sampled grid, then PCHIP-interpolate to `t`.
`solve(tblock)` must return the g-function evaluated on the uniform vector `tblock`.
# Arguments
    - `t`: Time vector (nt × 1) [s]
    - `solve`: Function that computes the g-function on a uniform time vector
# Output
    - `g`: g-function evaluated at `t` [°Cm/W]
"""
function _interp_gfunction(t, solve)
    t0, tend = float(minimum(t)), float(maximum(t))
    ts = Float64[]
    gs = Float64[]

    # Append the sub-sampled blocks that overlap the requested time range.
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

"""
    _apply_interp(t, interp::Bool, solve)

Optionally route a backend `solve(t)` through the sub-sampling grid. `interp = true` computes on the
reduced constant-step grid and interpolates to `t`; `interp = false` computes directly on `t`. Used
by `ground_response` to give one uniform keyword across every backend.
# Arguments
    - `t`: Time vector (nt × 1) [s]
    - `interp`: Whether to sub-sample and interpolate (default `true`)
    - `solve`: Function that computes the g-function on a uniform time vector
# Output
    - `g`: g-function evaluated at `t` [°Cm/W]
"""
function _apply_interp(t, interp::Bool, solve)
    return interp ? _interp_gfunction(t, solve) : solve(t)
end
