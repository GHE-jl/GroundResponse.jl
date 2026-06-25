# Exemple of use of the g_short_term function

using CairoMakie
include("../../src/g_short_term.jl")

# 1 - Input parameters
ks=2.13
Cs=2.0e6
kg=1.65
Cg=2.0e6
kp=0.4
Cp=1.9e6
Cf=4.2e6
ri=0.017
ro=0.022
rb=0.08
H=150
V̇=23.7/1000/60
D=0.029
dt=15 
tf = 3*24*3600 

# 2 - Transfer function

# validation of the input parameters
ks, Cs, kg, Cg, rb, H, V̇, D = input_validation(ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb, H, V̇, D, dt, tf)

# creating the Entering Water Temperature's timestep and the g function
t_EWT, g_EWT, _ = g_short_term(ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb, H, V̇, D, dt, tf)

# 3 - Creation of the figure

fig = Figure()
ax = Axis(
    fig[1, 1], 
    title = "g function computed with an Artificial Neural Network", 
    xlabel = "t/t_s (-)", 
    xscale = log10,
    ylabel = "g (-)" 
    )

lines!(
    ax, 
    t_EWT/(H^2/(9*ks/Cs)), 
    g_EWT, 
    color = :blue 
    )
    
fig