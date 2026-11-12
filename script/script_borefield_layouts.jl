# Visualization of all borefield layout configurations (~50 boreholes each).

import Pkg; Pkg.activate(@__DIR__)

using CairoMakie
using GroundResponse

B = 5.0    # Borehole spacing for all layouts [m]
R = 40.0   # Circle radius [m] — arc spacing ≈ 2π×40/50 ≈ 5.0 m

# Each entry: (subplot title, xy coordinates)
# Borehole counts are exact for all shapes.
layouts = [
    #  Shape                 Parameters            nb
    ("Rectangle  7×7 = 49",   borefield(:rectangle,      7,  7,  B)),  # 49
    ("Circle  nb = 50",        borefield(:circle,         50,    R)),   # 50
    ("Open rectangle  nb = 50",borefield(:open_rectangle, 15, 12, B)),  # 2*(15+12-2) = 50
    ("Line  nb = 50",          borefield(:line,           50,    B)),   # 50
    ("L-shape  26+25-1 = 50",  borefield(:L,              26, 25, B)),  # 50
    ("U-shape  18+2×16 = 50",  borefield(:U,              18, 17, B)),  # 18+2*16 = 50
]

f = Figure(size = (1050, 720))
Label(f[0, :], "Borefield layout configurations  (B = $(B) m)", fontsize = 16, font = :bold)

for (k, (title, xy)) in enumerate(layouts)
    row = (k - 1) ÷ 3 + 1
    col = mod1(k, 3)
    ax = Axis(f[row, col], title  = title, xlabel = "x (m)", ylabel = "y (m)",
        aspect = DataAspect())
    scatter!(ax, xy[:, 1], xy[:, 2], color = :steelblue, markersize = 8, strokecolor = :navy,
        strokewidth = 0.6)
end

display(f)

# Orientation of a borefield relative to the groundwater flow.
# The `angle` keyword rotates a layout counter-clockwise about the origin. Distances are unchanged,
# so the isotropic models (ILS/ICS/FLS) give the same g-function at any angle. The moving models
# (MILS/MFLS) do not, since the groundwater flow (+x) carries heat from upstream boreholes onto the
# downstream ones.

# Parameters for the example borefield and ground models
H, D, rb = 150.0, 4.0, 0.075 # Borehole depth, buried depth and radius [m]
ks, Cs = 2.5, 2.2e6         # Ground conductivity [W/mK] and volumetric heat capacity [J/m³K]
Cf, vD = 4.2e6, 1e-7        # Groundwater heat capacity [J/m³K] and Darcy velocity (~3 m/yr) [m/s]
B2 = 6.0                    # Spacing for this example [m]
m_flow  = MFLSModel(H, rb, D, ks, Cs, Cf, vD)
m_still = FLSModel(H, D, ks, Cs)
t = 3600.0 .* exp10.(range(log10(0.1), log10(8760 * 25), length = 300)) # 6 min → 25 years [s]
angles = 0:5:90

# Borefield layouts rotated by an angle α [°]
fields = [("Square 5×5", α -> borefield(:rectangle, 5, 5, B2; angle = α)),
          ("Rectangle 8×2", α -> borefield(:rectangle, 8, 2, B2; angle = α))]

# BC-II g-function of every layout at every angle, without (FLS) and with (MFLS) groundwater flow
sweep(f, m) = [ground_response(t, rb, f(α), m) for α in angles]
still = [sweep(f, m_still) for (_, f) in fields]
flow  = [sweep(f, m_flow) for (_, f) in fields]

# Figure: one row per layout; FLS (left) and MFLS (right) g-functions for every angle
f2 = Figure(size = (1000, 720))
Label(f2[0, 1:2], "Borefield orientation vs groundwater flow (+x, vD = $(vD) m/s)",
    fontsize = 16, font = :bold)
cmap, crange = :viridis, (first(angles), last(angles))
ty = t ./ (3600 * 8760)
for (k, (name, _)) in enumerate(fields)
    for (j, (gs, label)) in enumerate(((still[k], "FLS (no flow)"), (flow[k], "MFLS")))
        ax = Axis(f2[k, j], title = "$name — $label", xlabel = "Time (yr)",
            ylabel = "g-function (°Cm/W)", xscale = log10)
        for (α, g) in zip(angles, gs)
            lines!(ax, ty, g, color = α, colormap = cmap, colorrange = crange, linewidth = 1.5)
        end
    end
end
Colorbar(f2[1:2, 3], colormap = cmap, colorrange = crange, label = "Rotation angle (°)")

display(f2)
