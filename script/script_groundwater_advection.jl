# 2D spatial maps of the g-function at t = 10 yr for a single borehole at the origin.
# Compares FLS (radially symmetric, no advection) with MILS and MFLS (groundwater flow
# in the +x direction). The asymmetric plume elongation is clearly visible for MILS/MFLS.

import Pkg; Pkg.activate(@__DIR__)

using CairoMakie
using GroundResponse

# Parameters
H  = 150.0   # Borehole depth [m]
D  = 4.0     # Buried depth [m]
rb = 0.076   # Borehole radius [m]
ks = 3.0     # Ground thermal conductivity [W/mK]
Cs = 2.0e6   # Ground volumetric heat capacity [J/m³K]
Cf = 4.2e6   # Groundwater volumetric heat capacity [J/m³K]
vD = 1e-7    # Darcy velocity [m/s] — chosen for visible asymmetry at 10 yr

t_eval = 10.0 * 8760 * 3600.0   # 10 years [s]

# Spatial grid centered on the borehole
xs = range(-20.0, 40.0, step = 1.0)   # [m]  more extent downstream (+x)
ys = range(-25.0, 25.0, step = 1.0)   # [m]
nx_g, ny_g = length(xs), length(ys)

# Storage: g[i, j] = value at (x = xs[i], y = ys[j])
g_fls  = Matrix{Float64}(undef, nx_g, ny_g)
g_mils = Matrix{Float64}(undef, nx_g, ny_g)
g_mfls = Matrix{Float64}(undef, nx_g, ny_g)

# FLS — radially symmetric; vectorise over unique radii for speed
for (j, y) in enumerate(ys), (i, x) in enumerate(xs)
    r = max(sqrt(x^2 + y^2), rb)
    g_fls[i, j] = fls(t_eval, r, H, D, ks, Cs)
end

# MILS — analytical series, fast
for (j, y) in enumerate(ys), (i, x) in enumerate(xs)
    g_mils[i, j] = mils(t_eval, [x, y], rb, ks, Cs, Cf, vD)
end

# MFLS — numerical integration at each point, may take ~30 s
for (j, y) in enumerate(ys), (i, x) in enumerate(xs)
    g_mfls[i, j] = mfls(t_eval, [x, y], H, rb, D, ks, Cs, Cf, vD)
end

# Collect all values for consistent colour scaling across panels
all_vals = [g_fls; g_mils; g_mfls]

# Figure
f = Figure(size = (1150, 400))
Label(f[0, 1:3],
    "Groundwater advection effect — t = 10 yr, vD = $(vD) m/s, H = $(H) m",
    fontsize = 15, font = :bold)

panel_titles = [
    "FLS  (no advection)",
    "MILS  (infinite depth, vD = $(vD) m/s)",
    "MFLS  (finite depth, vD = $(vD) m/s)",
]
maps = [g_fls, g_mils, g_mfls]

for (k, (title, gmap)) in enumerate(zip(panel_titles, maps))
    ax = Axis(f[1, k],
        title  = title,
        xlabel = "x (m)",
        ylabel = k == 1 ? "y (m)" : "",
        aspect = DataAspect())

    hm = heatmap!(ax, collect(xs), collect(ys), gmap,
        colormap   = :turbo,
        colorrange = (minimum(vec(all_vals)), maximum(vec(all_vals))))

    # Borehole location marker
    scatter!(ax, [0.0], [0.0],
        marker      = :circle,
        color       = :white,
        strokecolor = :black,
        strokewidth = 1.5,
        markersize  = 9)

    # Groundwater flow direction arrow (downstream panels only)
    if k > 1
        arrows!(ax, [25.0], [-20.0], [8.0], [0.0],
            color     = :white,
            linewidth = 2,
            arrowsize = 10)
        text!(ax, 23.5, -17.5, text = "vD", color = :white, fontsize = 11)
    end

    # Shared colourbar on the right of the last panel
    k == 3 && Colorbar(f[1, 4], hm, label = "g-function (°Cm/W)")
end

display(f)
