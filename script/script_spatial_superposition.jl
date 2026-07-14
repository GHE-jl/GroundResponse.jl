# Spatial superposition: the boundary-condition hierarchy BC-I → BC-II → BC-III.
#
#   1. Boundary conditions — BC-I (uniform_flux), BC-II (bloc_matrix / successive_flux) and
#      BC-III (segment_response) on one 3×3 FLS borefield. The curves bracket the true response,
#      ordered BC-I ≥ BC-II ≥ BC-III.
#   2. BC-III convergence — the segment g-function as the number of segments per borehole grows.
#   3. Borefield size scaling — BC-II g-functions from 1×1 to 5×5.
#
# All three boundary conditions are available in the package and are selected through the optional
# `nseg` field of `FLSModel` (nseg = 1 → BC-I/BC-II, nseg > 1 → BC-III).

import Pkg; Pkg.activate(@__DIR__)

using CairoMakie
using GroundResponse

# Common ground parameters
H  = 150.0    # Borehole depth [m]
D  = 4.0      # Buried depth [m]
rb = 0.076    # Borehole radius [m]
ks = 3.0      # Ground thermal conductivity [W/mK]
Cs = 2.0e6    # Ground volumetric heat capacity [J/m³K]
B  = 5.0      # Borehole spacing [m]

# Log-spaced time vector: 1 h → 25 yr. Kept moderate (nt = 60) because the BC-III block solve grows
# with (nb · nseg · nt); this resolution is plenty for a smooth g-function.
t  = 3600.0 .* exp10.(range(0, log10(8760 * 25), length = 60))
t̃ = t ./ (3600 * 8760)   # Time in years

m_fls = FLSModel(H, D, ks, Cs)
xy33  = borefield(:rectangle, 3, 3, B)

# 1. Boundary-condition hierarchy on a 3×3 FLS borefield
println("Computing BC-I / BC-II / BC-III g-functions (3×3 FLS)...")
g_bc1 = uniform_flux(t, rb, xy33, m_fls)                       # BC-I  — uniform flux
g_bc2 = bloc_matrix(t, rb, xy33, m_fls)                        # BC-II — equal mean temperature
g_sf  = successive_flux(t, rb, xy33, m_fls)                    # BC-II — iterative (same BC as bloc)
g_bc3 = segment_response(t, rb, xy33, FLSModel(H, D, ks, Cs, 8))  # BC-III — uniform wall temperature

# 2. BC-III convergence with the number of segments
println("Computing BC-III convergence with nseg...")
nsegs = [1, 2, 4, 8, 16]
gs_nseg = [segment_response(t, rb, xy33, FLSModel(H, D, ks, Cs, ns)) for ns in nsegs]

# 3. Borefield size scaling (BC-II — fast for larger fields)
println("Computing BC-II g-functions for borefield size scaling...")
sizes = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]
labels_size = ["1×1 = 1", "2×2 = 4", "3×3 = 9", "4×4 = 16", "5×5 = 25"]
gs_size = [ground_response(t, rb, borefield(:rectangle, nx, ny, B), m_fls) for (nx, ny) in sizes]

# Figures
# Figure 1: boundary-condition hierarchy (left) and BC-III convergence (right)
f1 = Figure(size = (1100, 480))

ax1 = Axis(f1[1, 1], xlabel = "Time (yr)", ylabel = "g-function (°Cm/W)",
    title = "Boundary conditions BC-I → BC-III  (3×3 FLS, B = $(B) m)", xscale = log10)
lines!(ax1, t̃, g_bc1, linewidth = 3, color = :darkgreen, label = "BC-I  (uniform flux)")
lines!(ax1, t̃, g_bc2, linewidth = 3, color = :steelblue, label = "BC-II  (bloc matrix)")
lines!(ax1, t̃, g_sf, linewidth = 2, linestyle = :dash, color = :royalblue,
    label = "BC-II  (successive flux)")
lines!(ax1, t̃, g_bc3, linewidth = 3, color = :tomato, label = "BC-III  (segment, nseg = 8)")
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
for (k, (g, lbl)) in enumerate(zip(gs_size, labels_size))
    lines!(ax3, t̃, g, linewidth = 2.5, color = Cycled(k), label = lbl)
end
axislegend(ax3, position = :lt)

display(f2)
