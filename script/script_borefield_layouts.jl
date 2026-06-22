# Visualization of all borefield layout configurations (~50 boreholes each).

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
    ax = Axis(f[row, col],
        title  = title,
        xlabel = "x (m)",
        ylabel = "y (m)",
        aspect = DataAspect())
    scatter!(ax, xy[:, 1], xy[:, 2],
        color        = :steelblue,
        markersize   = 8,
        strokecolor  = :navy,
        strokewidth  = 0.6)
end

display(f)
