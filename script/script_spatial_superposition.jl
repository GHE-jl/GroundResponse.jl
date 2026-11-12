# Spatial superposition: the boundary conditions BC-I → BC-II → BC-III, and the marching solver
# behind BC-II / BC-III compared with the legacy solvers.
#
#   1. Boundary conditions — bc_I, bc_II and bc_III on one 3×3 FLS borefield. The curves bracket the
#      true response, ordered BC-I ≥ BC-II ≥ BC-III.
#   2. BC-III convergence — the segment g-function as the number of segments per borehole grows.
#   3. Borefield size scaling — BC-II g-functions from 1×1 to 5×5.
#   4. Solvers — computing time of bc_II / bc_III against the legacy solvers on the same problem, and
#      their difference to it (the block matrix solves the same system at once; the successive flux
#      iterates towards it).
#   5. Time sampling — the default interp = true against a direct solve on an hourly grid.

import Pkg; Pkg.activate(@__DIR__)

using CairoMakie
using GroundResponse
using Printf

# Common ground parameters
H  = 150.0    # Borehole depth [m]
D  = 4.0      # Buried depth [m]
rb = 0.076    # Borehole radius [m]
ks = 3.0      # Ground thermal conductivity [W/mK]
Cs = 2.0e6    # Ground volumetric heat capacity [J/m³K]
B  = 5.0      # Borehole spacing [m]

# Log-spaced time vector: 1 h → 25 yr (solved on the internal sub-sampling grid, interp = true)
t  = 3600.0 .* exp10.(range(0, log10(8760 * 25), length = 60))
t̃ = t ./ (3600 * 8760)   # Time in years

m_fls = FLSModel(H, D, ks, Cs)
xy33  = borefield(:rectangle, 3, 3, B)

# 1. Boundary conditions on a 3×3 FLS borefield
g_bc1 = bc_I(t, rb, xy33, m_fls)                          # BC-I   — uniform flux
g_bc2 = bc_II(t, rb, xy33, m_fls)                         # BC-II  — equal mean temperature
g_bc3 = bc_III(t, rb, xy33, FLSModel(H, D, ks, Cs, 8))    # BC-III — uniform wall temperature

# 2. BC-III convergence with the number of segments
nsegs   = [1, 2, 4, 8, 16]
gs_nseg = [bc_III(t, rb, xy33, FLSModel(H, D, ks, Cs, ns)) for ns in nsegs]

# 3. Borefield size scaling (BC-II)
sizes   = [1, 2, 3, 4, 5]
gs_size = [ground_response(t, rb, borefield(:rectangle, n, n, B), m_fls) for n in sizes]

# 4. Solvers: computing time and difference to bc_II / bc_III on the same problem
timed(f) = (f(); minimum(@elapsed(f()) for _ in 1:2))     # first call compiles
reldiff(a, b) = maximum(abs.(a .- b) ./ b)

println("BC-II  (FLS, n×n field, t = 1 h → 25 yr)")
@printf("%5s | %9s | %22s | %22s\n", "nb", "bc_II [s]", "successive_legacy", "bloc_matrix_legacy")
bc2_runs = []
for n in (2, 3, 5, 8, 10, 15, 20)
    xy = borefield(:rectangle, n, n, B)
    ref = bc_II(t, rb, xy, m_fls)
    run = (nb = n^2, new = timed(() -> bc_II(t, rb, xy, m_fls)),
        sf = timed(() -> successive_flux_legacy(t, rb, xy, m_fls)),
        dsf = reldiff(successive_flux_legacy(t, rb, xy, m_fls), ref),
        bm = n ≤ 8 ? timed(() -> bloc_matrix_legacy(t, rb, xy, m_fls)) : NaN,
        dbm = n ≤ 8 ? reldiff(bloc_matrix_legacy(t, rb, xy, m_fls), ref) : NaN)
    push!(bc2_runs, run)
    @printf("%5d | %9.3f | %8.3f s  (Δ %.1e) | %8.3f s  (Δ %.1e)\n", run.nb, run.new, run.sf,
        run.dsf, run.bm, run.dbm)
end

println("\nBC-III (FLS, n×n field, 8 segments per borehole)")
@printf("%5s | %10s | %22s | %22s\n", "nS", "bc_III [s]", "marching_legacy", "segment_response_legacy")
bc3_runs = []
for n in (2, 3, 5, 8)
    xy = borefield(:rectangle, n, n, B); m8 = FLSModel(H, D, ks, Cs, 8)
    ref = bc_III(t, rb, xy, m8)
    run = (nS = 8n^2, new = timed(() -> bc_III(t, rb, xy, m8)),
        mr = timed(() -> segment_response_marching_legacy(t, rb, xy, m8)),
        dmr = reldiff(segment_response_marching_legacy(t, rb, xy, m8), ref),
        bm = n ≤ 3 ? timed(() -> segment_response_legacy(t, rb, xy, m8)) : NaN,
        dbm = n ≤ 3 ? reldiff(segment_response_legacy(t, rb, xy, m8), ref) : NaN)
    push!(bc3_runs, run)
    @printf("%5d | %10.3f | %8.3f s  (Δ %.1e) | %8.3f s  (Δ %.1e)\n", run.nS, run.new, run.mr,
        run.dmr, run.bm, run.dbm)
end

# 5. Time sampling: interp = true (sub-sampled blocks, PCHIP in log t) against a direct solve on an
# hourly grid over one year (interp = false needs t = Δt, 2Δt, …)
th = 3600.0 .* (1:8760)
gh = bc_II(th, rb, xy33, m_fls; interp = false)
gi = bc_II(th, rb, xy33, m_fls)
@printf("\nTime sampling, 3×3 BC-II over 1 yr: interp = true vs hourly direct solve, max Δ %.1e\n",
    reldiff(gi, gh))

# Figures
# Figure 1: boundary conditions (left) and BC-III convergence (right)
f1 = Figure(size = (1100, 480))
ax1 = Axis(f1[1, 1], xlabel = "Time (yr)", ylabel = "g-function (°Cm/W)",
    title = "Boundary conditions BC-I → BC-III  (3×3 FLS, B = $(B) m)", xscale = log10)
lines!(ax1, t̃, g_bc1, linewidth = 3, color = :darkgreen, label = "BC-I  (uniform flux)")
lines!(ax1, t̃, g_bc2, linewidth = 3, color = :steelblue, label = "BC-II  (equal mean temperature)")
lines!(ax1, t̃, g_bc3, linewidth = 3, color = :tomato, label = "BC-III  (nseg = 8)")
axislegend(ax1, position = :lt)

ax2 = Axis(f1[1, 2], xlabel = "Time (yr)", ylabel = "g-function (°Cm/W)",
    title = "BC-III convergence with segments per borehole (3×3)", xscale = log10)
for (k, (g, ns)) in enumerate(zip(gs_nseg, nsegs))
    lines!(ax2, t̃, g, linewidth = 2.5, color = Cycled(k),
        label = ns == 1 ? "nseg = 1 (= BC-II)" : "nseg = $ns")
end
axislegend(ax2, position = :lt)
display(f1)

# Figure 2: borefield size scaling
f2 = Figure(size = (620, 480))
ax3 = Axis(f2[1, 1], xlabel = "Time (yr)", ylabel = "g-function (°Cm/W)",
    title = "Borefield size scaling — FLS + BC-II, B = $(B) m", xscale = log10)
for (k, (g, n)) in enumerate(zip(gs_size, sizes))
    lines!(ax3, t̃, g, linewidth = 2.5, color = Cycled(k), label = "$n×$n = $(n^2)")
end
axislegend(ax3, position = :lt)
display(f2)

# Figure 3: computing time of the solvers
f3 = Figure(size = (1100, 460))
ax4 = Axis(f3[1, 1], xlabel = "Boreholes", ylabel = "Computing time (s)", xscale = log10,
    yscale = log10, title = "BC-II")
nb = [r.nb for r in bc2_runs]
scatterlines!(ax4, nb, [r.new for r in bc2_runs], linewidth = 2.5, label = "bc_II (marching)")
scatterlines!(ax4, nb, [r.sf for r in bc2_runs], linewidth = 2, label = "successive_flux_legacy")
ok = isfinite.([r.bm for r in bc2_runs])
scatterlines!(ax4, nb[ok], [r.bm for r in bc2_runs][ok], linewidth = 2, label = "bloc_matrix_legacy")
axislegend(ax4, position = :lt)

ax5 = Axis(f3[1, 2], xlabel = "Segments (boreholes × 8)", ylabel = "Computing time (s)",
    xscale = log10, yscale = log10, title = "BC-III")
nS = [r.nS for r in bc3_runs]
scatterlines!(ax5, nS, [r.new for r in bc3_runs], linewidth = 2.5, label = "bc_III (marching)")
scatterlines!(ax5, nS, [r.mr for r in bc3_runs], linewidth = 2,
    label = "segment_response_marching_legacy")
ok = isfinite.([r.bm for r in bc3_runs])
scatterlines!(ax5, nS[ok], [r.bm for r in bc3_runs][ok], linewidth = 2,
    label = "segment_response_legacy")
axislegend(ax5, position = :lt)
display(f3)
