"""
Script testing the nodes generation to apply to analytical models.
"""

using BenchmarkTools
using CairoMakie

includet("../src/GroundHeatExchanger.jl")
using .GroundHeatExchanger

# Define paremeters
ti, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()
t = 1:3600.0:ti[end]        # [s]

# Generate nodes and evaluate the ILS model at those nodes
s = set_nodes(length(t), 100)
@time ils(t[s], kg, Cg, rb)
gᵢ = ils(t[s], kg, Cg, rb)

# Evaluate the ILS model at the original time vector
@time ils(t, kg, Cg, rb)
g = ils(t, kg, Cg, rb)

# Interpolate the ILS model at the original time vector using PCHIP interpolation
@time pchip_interpolation(t[s], gᵢ, t)
g̃ = pchip_interpolation(t[s], gᵢ, t)

# Compute RMSE between the original ILS model and the interpolated version
rmse = sqrt(mean((g - g̃).^2))

# Figure
fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time (s)", ylabel="ILS-model (°Cm/W)", xscale=log10)
lines!(ax, t, g, label="ILS")
lines!(ax, t, g̃, label="PCHIP")
scatter!(ax, t[s], gᵢ, markersize=3, label="Nodes")
text!(ax, 0.1, 0.5; text="RMSE = $(round(rmse, digits=3))", space=:relative, align=(:left, :center))
axislegend(ax, position=:lt)
display(fig)