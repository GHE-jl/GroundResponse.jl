# Single-borehole g-function comparison: ILS, ICS, FLS, MILS, MFLS.

import Pkg; Pkg.activate(@__DIR__)

using CairoMakie
using GroundResponse

# Ground and borehole parameters
H  = 150.0   # Borehole depth [m]
D  = 4.0     # Buried depth [m]
rb = 0.076   # Borehole radius [m]
ks = 3.0     # Ground thermal conductivity [W/mK]
Cs = 2.0e6   # Ground volumetric heat capacity [J/m³K]
Cf = 4.2e6   # Groundwater volumetric heat capacity [J/m³K]
vD = 1e-7    # Darcy velocity [m/s]

# Log-spaced time vector: 0.1 h → 25 yr
t = 3600.0 .* exp10.(range(log10(0.1), log10(8760 * 25), length = 300))

# Single-borehole g-functions at the borehole wall
g_ils  = ils(t, rb, ks, Cs)
g_ics  = ics(t, rb, rb, ks, Cs)
g_fls  = fls(t, rb, H, D, ks, Cs)
g_mils = mils(t, [0.0, 0.0], rb, ks, Cs, Cf, vD)
g_mfls = mfls(t, [0.0, 0.0], H, rb, D, ks, Cs, Cf, vD)

t̃ = t ./ (3600 * 8760)   # Time in years

f = Figure(size = (720, 460))
ax = Axis(f[1, 1],
    xlabel = "Time (yr)",
    ylabel = "g-function (°Cm/W)",
    title  = "Single borehole — all models  (H = $(H) m, rb = $(rb) m)",
    xscale = log10)
lines!(ax, t̃, g_ils,  linewidth = 3, linestyle = :solid,       label = "ILS")
lines!(ax, t̃, g_ics,  linewidth = 3, linestyle = :dash,        label = "ICS")
lines!(ax, t̃, g_fls,  linewidth = 3, linestyle = :dot,         label = "FLS")
lines!(ax, t̃, g_mils, linewidth = 3, linestyle = :dashdot,     label = "MILS  (vD = $(vD) m/s)")
lines!(ax, t̃, g_mfls, linewidth = 3, linestyle = :dashdotdot,  label = "MFLS  (vD = $(vD) m/s)")
axislegend(ax, position = :lt)
display(f)
