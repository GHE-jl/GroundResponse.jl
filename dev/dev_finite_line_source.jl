"""
Script that tests the multiple implementation of the finite line source model to check if they
provide the same results.
"""

using BenchmarkTools
using CairoMakie

includet("../src/ground_models/finite_line_source.jl")
includet("../src/Utils.jl")

# Define paremeters
t, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()
nt = length(t)

# Borefield definition
xy = borefield_xy(2, 2, 5.)
r, rᵥ, rᵤ, rᵢ, θ, nb = borefield_radius(xy, rb)

# Finite line source for multiple time steps and 1 radius
fls1 = fls(t[end], ks, Cs, rb, H, D)
fls2 = fls(t, ks, Cs, rb, H, D)
fls2[end] == fls1

# Finite line source for 1 time step and multiple radius
fls3 = fls(t[end], ks, Cs, rᵤ[end], H, D)
fls4 = fls(t[end], ks, Cs, rᵤ, H, D)
fls4[end] == fls3

# Finite line source for multiple time steps and multiple radius
fls5 = fls(t, ks, Cs, rᵤ, H, D)
fls6 = fls(t, ks, Cs, r, H, D)