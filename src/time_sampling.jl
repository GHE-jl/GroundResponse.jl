# The constant-step sub-sampling grid and the PCHIP interpolation used to (a) satisfy the
# constant-step requirement of the BC-II / BC-III solvers and (b) bound the cost of any backend on a
# large `t`.

# Fixed sub-sampling scheme for the spatial-superposition solvers of BC-II and BC-III (`bc_II`,
# `bc_III` and the legacy solvers). Their load history reads the response at the lag
# `t_p − t_{l−1}` as `g[p − l + 1]`, which holds only on a uniform grid `t = Δt, 2Δt, …, nΔt`, so an
# arbitrary (e.g. log-spaced) `t` would silently corrupt the temporal superposition. The g-function
# is instead solved on the blocks below, each uniform from t = 0, following the geometric
# hour→decade progression of Cimmino et al. (2013) and Nguyen & Pasquier (2021). Each block restarts
# from t = 0 with its own step, so the flux redistribution during its first step is not resolved;
# the quarterly and yearly blocks keep that error small at the time scales where the boreholes
# interact (about 0.2 % for 100 boreholes, against 0.9 % with yearly steps to 10 years and decades
# beyond). The nodes are then interpolated by PCHIP in log t, the natural variable of a g-function.
# The same routine also serves BC-I and the single-borehole models, where it is a pure performance
# approximation rather than a correctness requirement.

const _SUBSAMPLE_BLOCKS = ((2.0,          30),   # every 2 s   → 60 s
                           (120.0,        30),   # every 2 min → 1 h
                           (3600.0,       24),   # every 1 h   → 1 day
                           (86_400.0,     30),   # every 1 day → 1 month (30-day)
                           (2_592_000.0,  12),   # every 1 month → 1 year (360-day)
                           (7_776_000.0,  40),   # every 3 months → 10 years
                           (31_104_000.0, 50))   # every 1 year → 50 years

const _SUBSAMPLE_END    = 1_555_200_000.0        # end of the last block (50 years) [s]
const _SUBSAMPLE_DECADE = 311_040_000.0          # 1 decade [s]; beyond 50 years, decades as needed

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

Solve a temporal g-function on the fixed sub-sampled grid, then PCHIP-interpolate to `t` in log t.
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
    if tend > _SUBSAMPLE_END               # beyond 50 years → add decade-spaced points
        ndec = ceil(Int, tend / _SUBSAMPLE_DECADE)
        tb = collect(range(_SUBSAMPLE_DECADE, ndec * _SUBSAMPLE_DECADE; length = ndec))
        _append_strictly_increasing!(ts, gs, tb, solve(tb))
    end
    @assert length(ts) ≥ 2 "sub-sampling produced too few nodes; pass interp = false"
    itp = Interpolator(log.(ts), gs)
    return [itp(log(clamp(tt, ts[1], ts[end]))) for tt in t]
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

"""
    _check_uniform_steps(t)

Reject a `t` that is not `Δt, 2Δt, …, nΔt` when a BC-II / BC-III solver is asked to solve directly
on it (`interp = false`): its load history would be superposed at the wrong lags.
# Arguments
    - `t`: Time value or vector [s]
"""
function _check_uniform_steps(t)
    Δt = first(t)
    all(k -> isapprox(t[k], k * Δt; rtol = 1e-9), eachindex(t)) || throw(ArgumentError(
        "interp = false requires uniform time steps from t = 0 (t = Δt, 2Δt, …, nΔt); " *
        "use interp = true for any other t"))
    return nothing
end
