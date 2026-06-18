"""
Script testing the spatial superposition with the finite line source (FLS) model.
"""

using BenchmarkTools
using CairoMakie

includet("../src/GroundHeatExchanger.jl")
using .GroundHeatExchanger

# Define paremeters
t, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()

# Define the borefield geometry
nx, ny, B = 5, 5, 5.
xy = borefield_xy(nx, ny, B)
r, rᵤ, rᵥ, rᵢ, θ, nb = borefield_radius(xy, rb)

# Define 3D matrix of g-function for different radius of the borefield
g = fls(t, H, r, D, ks, Cs)

# Spatial superposition with bloc matrix (around 100x slower than successive flux)
@time g_bm1 = bloc_matrix(g)
@time g_bm2 = bloc_matrix(t, H, rb, D, ks, Cs, xy)

# Spatial superposition with successive flux
@time g_sf1 = successive_flux(g)
@time g_sf2 = successive_flux(t, H, rb, D, ks, Cs, xy)

# Figure
fig = Figure()
ax = Axis(fig[1, 1], xlabel = L"$t$ (s)", ylabel = L"$g$ (°Cm/W)", xscale = log10)
lines!(ax, t, g_bm1, color = "red", linewidth=3, linestyle = :solid, label = "Bloc matrix $nx x $ny")
lines!(ax, t, g_bm2, color = "purple", linewidth=2.5, linestyle = :dashdot, label = "Bloc matrix $nx x $ny")
lines!(ax, t, g_sf1, color = "blue", linewidth=2, linestyle = :dash, label = "Successive flux $nx x $ny")
lines!(ax, t, g_sf2, color = "green", linewidth=1.5, linestyle = :dot, label = "Successive flux $nx x $ny")
axislegend(ax, position = :lt)
display(fig)