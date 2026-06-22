# Spatial superposition: bloc matrix vs successive flux for FLS borefields.
#   1. Method agreement — compare both algorithms on a 3×3 borefield.
#   2. Borefield size scaling — g-functions from 1×1 to 5×5.

using CairoMakie
using GroundResponse

# Parameters
H  = 150.0   # Borehole depth [m]
D  = 4.0     # Buried depth [m]
rb = 0.076   # Borehole radius [m]
ks = 3.0     # Ground thermal conductivity [W/mK]
Cs = 2.0e6   # Ground volumetric heat capacity [J/m³K]
B  = 5.0     # Borehole spacing [m]

# Log-spaced time vector: 1 h → 25 yr
t = 3600.0 .* exp10.(range(0, log10(8760 * 25), length = 200))
t̃ = t ./ (3600 * 8760)   # Time in years

m = FLSModel(H, D, ks, Cs)

# 1. Method comparison on a 5×5 borefield
xy55  = borefield(:rectangle, 5, 5, B)

g_bm  = bloc_matrix(t, rb, xy55, m)
g_sf  = successive_flux(t, rb, xy55, m)

rel_err = maximum(abs, (g_bm .- g_sf) ./ max.(abs.(g_bm), 1e-12))
println("Max relative error (bloc matrix vs successive flux): ",
        round(rel_err * 100, digits = 4), " %")

# 2. Borefield size scaling (successive_flux via ground_response)
sizes  = [(1,1), (2,2), (3,3), (4,4), (5,5)]
labels = ["1 borehole", "2×2 = 4", "3×3 = 9", "4×4 = 16", "5×5 = 25"]
println("\nComputing g-functions for borefield size scaling...")
gs = [ground_response(t, rb, borefield(:rectangle, nx, ny, B), m) for (nx, ny) in sizes]

# Figure
f = Figure(size = (1100, 480))

# Left: method agreement
ax1 = Axis(f[1, 1],
    xlabel = "Time (yr)", ylabel = "g-function (°Cm/W)",
    title  = "Bloc matrix vs successive flux  (5×5 FLS, B = $(B) m)",
    xscale = log10)
lines!(ax1, t̃, g_bm, linewidth = 3, color = :steelblue,
    label = "Bloc matrix")
lines!(ax1, t̃, g_sf, linewidth = 2, color = :tomato, linestyle = :dash,
    label = "Successive flux")
axislegend(ax1, position = :lt)

# Right: size scaling
ax2 = Axis(f[1, 2],
    xlabel = "Time (yr)", ylabel = "g-function (°Cm/W)",
    title  = "Borefield size scaling — FLS, B = $(B) m",
    xscale = log10)
colors = Makie.wong_colors()
for (k, (g, lbl)) in enumerate(zip(gs, labels))
    lines!(ax2, t̃, g, linewidth = 2.5, color = colors[k], label = lbl)
end
axislegend(ax2, position = :lt)

display(f)
