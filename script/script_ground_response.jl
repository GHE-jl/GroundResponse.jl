# High-level ground_response interface: AbstractGroundModel dispatch, benchmark across
# models (single borehole) and borefield sizes (FLS, successive flux).

import Pkg; Pkg.activate(@__DIR__)

using BenchmarkTools
using CairoMakie
using GroundResponse

# Parameters
H  = 150.0   # Borehole depth [m]
D  = 4.0     # Buried depth [m]
rb = 0.076   # Borehole radius [m]
ks = 3.0     # Ground thermal conductivity [W/mK]
Cs = 2.0e6   # Ground volumetric heat capacity [J/m³K]
Cf = 4.2e6   # Groundwater volumetric heat capacity [J/m³K]
vD = 1e-7    # Darcy velocity [m/s]
B  = 5.0     # Borehole spacing [m]

# 300 log-spaced time steps: 1 h → 25 yr
t = 3600.0 .* exp10.(range(0, log10(8760 * 25), length = 300))

# Ground models — all subtypes of AbstractGroundModel
m_ils  = ILSModel(ks, Cs)
m_ics  = ICSModel(rb, ks, Cs)
m_fls  = FLSModel(H, D, ks, Cs)
m_mils = MILSModel(rb, ks, Cs, Cf, vD)
m_mfls = MFLSModel(H, rb, D, ks, Cs, Cf, vD)

# Single borehole: 1×2 coordinate matrix with the borehole at the origin
xy1 = [0.0 0.0]

# Single borehole benchmark
println("  Single borehole  —  ground_response  (nt = $(length(t)))")
print("  ILS  : "); @btime ground_response(t, rb, xy1, m_ils)
print("  ICS  : "); @btime ground_response(t, rb, xy1, m_ics)
print("  FLS  : "); @btime ground_response(t, rb, xy1, m_fls)
print("  MILS : "); @btime ground_response(t, rb, xy1, m_mils)
print("  MFLS : "); @btime ground_response(t, rb, xy1, m_mfls)

# Borefield size scaling (FLS, successive flux)
println("  FLS borefield  —  ground_response  (successive flux)")
bfield_sizes = [(2, 2), (3, 3), (4, 4), (5, 5)]
for (nx, ny) in bfield_sizes
    xy = borefield(:rectangle, nx, ny, B)
    nb = size(xy, 1)
    print("  $(nx)×$(ny) (nb = $(nb)): ")
    @btime ground_response($t, $rb, $xy, $m_fls)
end

# Figure: g-functions for different borefield sizes
t̃ = t ./ (3600 * 8760)   # Time in years

f = Figure(size = (800, 500))
ax = Axis(f[1, 1],
    xlabel = "Time (yr)", ylabel = "g-function (°Cm/W)",
    title  = "FLS borefield g-functions via ground_response  (B = $(B) m)",
    xscale = log10)

all_sizes  = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]
all_labels = ["Single borehole", "2×2 = 4", "3×3 = 9", "4×4 = 16", "5×5 = 25"]
colors     = Makie.wong_colors()

for (k, ((nx, ny), lbl)) in enumerate(zip(all_sizes, all_labels))
    xy = borefield(:rectangle, nx, ny, B)
    g  = ground_response(t, rb, xy, m_fls)
    lines!(ax, t̃, g, linewidth = 2.5, color = colors[k], label = lbl)
end
axislegend(ax, position = :lt)
display(f)
