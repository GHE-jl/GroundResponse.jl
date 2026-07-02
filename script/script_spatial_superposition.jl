# Spatial superposition: method agreement, borefield size scaling, and model comparison.
#   1. Method agreement — compare bloc_matrix vs successive_flux on a 5×5 FLS borefield.
#   2. Borefield size scaling — g-functions from 1×1 to 5×5 (FLS).
#   3. Model comparison — ILS, ICS, FLS, MILS, MFLS on a 3×3 borefield.

import Pkg; Pkg.activate(@__DIR__)

using CairoMakie
using GroundResponse

# Common ground parameters
H  = 150.0   # Borehole depth [m]
D  = 4.0     # Buried depth [m]
rb = 0.076   # Borehole radius [m]
ks = 3.0     # Ground thermal conductivity [W/mK]
Cs = 2.0e6   # Ground volumetric heat capacity [J/m³K]
Cf = 4.18e6  # Groundwater volumetric heat capacity [J/m³K]
vD = 1e-6    # Darcy velocity [m/s]
B  = 5.0     # Borehole spacing [m]

# Log-spaced time vector: 1 h → 25 yr
t  = 3600.0 .* exp10.(range(0, log10(8760 * 25), length = 200))
t̃  = t ./ (3600 * 8760)   # Time in years

# 1. Method agreement on a 5×5 FLS borefield
m_fls = FLSModel(H, D, ks, Cs)
xy55  = borefield(:rectangle, 5, 5, B)

g_bm  = bloc_matrix(t, rb, xy55, m_fls)
g_sf  = successive_flux(t, rb, xy55, m_fls)

rel_err = maximum(abs, (g_bm .- g_sf) ./ max.(abs.(g_bm), 1e-12))
println("Max relative error (bloc matrix vs successive flux): ",
        round(rel_err * 100, digits = 4), " %")

# 2. Borefield size scaling
sizes       = [(1,1), (2,2), (3,3), (4,4), (5,5)]
labels_size = ["1×1 = 1", "2×2 = 4", "3×3 = 9", "4×4 = 16", "5×5 = 25"]
println("\nComputing g-functions for borefield size scaling...")
gs_size = [ground_response(t, rb, borefield(:rectangle, nx, ny, B), m_fls)
           for (nx, ny) in sizes]

# 3. Model comparison on a 3×3 borefield
println("\nComputing g-functions for model comparison (3×3 borefield)...")
xy33 = borefield(:rectangle, 3, 3, B)
models = [
    ILSModel(ks, Cs),
    ICSModel(rb, ks, Cs),
    FLSModel(H, D, ks, Cs),
    # MILSModel(rb, ks, Cs, Cf, vD),
    MFLSModel(H, rb, D, ks, Cs, Cf, vD),
]
# model_labels = ["ILS", "ICS", "FLS", "MILS", "MFLS"]#TODO: Fix for the MILS
model_labels = ["ILS", "ICS", "FLS", "MFLS"]
gs_models = [successive_flux(t, rb, xy33, m) for m in models]

# Figures
# Figure 1: method agreement (left) and size scaling (right)
f1 = Figure(size = (1100, 480))

ax1 = Axis(f1[1, 1],
    xlabel = "Time (yr)", ylabel = "g-function (°Cm/W)",
    title  = "Bloc matrix vs successive flux  (5×5 FLS, B = $(B) m)",
    xscale = log10)
lines!(ax1, t̃, g_bm, linewidth = 3, color = :steelblue, label = "Bloc matrix")
lines!(ax1, t̃, g_sf, linewidth = 2, color = :tomato, linestyle = :dash,
       label = "Successive flux")
axislegend(ax1, position = :lt)

ax2 = Axis(f1[1, 2],
    xlabel = "Time (yr)", ylabel = "g-function (°Cm/W)",
    title  = "Borefield size scaling — FLS, B = $(B) m",
    xscale = log10)
for (k, (g, lbl)) in enumerate(zip(gs_size, labels_size))
    lines!(ax2, t̃, g, linewidth = 2.5, color = colors[k], label = lbl)
end
axislegend(ax2, position = :lt)

display(f1)

# Figure 2: model comparison on a 3×3 borefield
f2 = Figure(size = (800, 500))

ax3 = Axis(f2[1, 1],
    xlabel = "Time (yr)", ylabel = "g-function (°Cm/W)",
    title  = "Ground model comparison — 3×3 borefield  " *
             "(B = $(B) m, H = $(H) m, D = $(D) m, ks = $(ks) W/mK)\n" *
             "Moving models: vD = $(vD) m/s, Cf = $(Cf/1e6) MJ/m³K",
    xscale = log10)
linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
for (k, (g, lbl, ls)) in enumerate(zip(gs_models, model_labels, linestyles))
    lines!(ax3, t̃, g, linewidth = 2.5, color = colors[k], linestyle = ls, label = lbl)
end
axislegend(ax3, position = :lt)

display(f2)
