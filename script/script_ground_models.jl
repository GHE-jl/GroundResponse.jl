"""
Script showcasing the analytical models from GHEModels.jl using a single borehole.
"""

using BenchmarkTools
using CairoMakie

includet("../src/GroundHeatExchanger.jl")
using .GroundHeatExchanger

# Define paremeters
t, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()

# Run models
@time g_ils = ils(t, ks, Cs, rb)
@time g_ics = ics(t, ks, Cs, rb)
@time g_fls = fls(t, ks, Cs, rb, H, D)
@time g_mils = mils(t, ks, Cs, Cf, rb, vD)
@time g_mfls = mfls(t, ks, Cs, Cf, [0, 0], rb, H, D, vD)

t̃ = t / (3600 * 24 * 365)
f = Figure()
ax = Axis(f[1, 1], xlabel = "Time (yr)", ylabel = "g-function (°Cm/W)", xscale = log10)
lines!(ax, t̃, g_ils, linewidth = 5, linestyle = :solid, label = "ILS")
lines!(ax, t̃, g_ics, linewidth = 4, linestyle = :dash, label = "ICS")
lines!(ax, t̃, g_fls, linewidth = 3, linestyle = :dot, label = "FLS")
lines!(ax, t̃, g_mils, linewidth = 2, linestyle = :dashdot, label = "MILS")
lines!(ax, t̃, g_mfls, linewidth = 1, linestyle = :dashdotdot, label = "MFLS")
axislegend(ax, position = :lt)
display(f);