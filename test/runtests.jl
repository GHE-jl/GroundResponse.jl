using Test
using GroundResponse
using LinearAlgebra: diag

# Shared physical parameters used across all test sets
const rb  = 0.076    # borehole / cylinder radius [m]
const H   = 150.0    # borehole depth [m]
const D   = 4.0      # buried depth [m]
const ks  = 3.0      # ground thermal conductivity [W/mK]
const Cs  = 2.0e6    # ground volumetric heat capacity [J/m³K]
const Cf  = 4.2e6    # groundwater volumetric heat capacity [J/m³K]
const vD  = 1e-6     # Darcy velocity [m/s]
const B   = 5.0      # borehole spacing [m]

# ILS — Infinite Line Source
@testset "ILS — Infinite Line Source" begin

    @testset "Scalar time and radius" begin
        g = ils(3600.0, rb, ks, Cs)
        @test g > 0
        @test isfinite(g)
        @test isa(g, AbstractFloat)
    end

    @testset "Integer time is promoted to Float" begin
        @test isa(ils(3600, rb, ks, Cs), AbstractFloat)
    end

    @testset "Vector time → vector, monotonically increasing" begin
        t = 60.0:60.0:3600.0
        g = ils(t, rb, ks, Cs)
        @test isa(g, Vector)
        @test length(g) == length(t)
        @test all(isfinite, g)
        @test all(g .> 0)
        @test all(diff(g) .> 0)
    end

    @testset "Vector radius → vector, monotonically decreasing" begin
        r_vec = [0.076, 0.2, 0.5, 1.0]
        g = ils(3600.0, r_vec, ks, Cs)
        @test isa(g, Vector)
        @test length(g) == length(r_vec)
        @test all(isfinite, g)
        @test all(diff(g) .< 0)
    end

    @testset "Vector time × vector radius → matrix" begin
        t = 3600.0 .* exp10.(range(0, log10(8760), length=8))
        r_vec = [0.076, 0.2, 0.5]
        g = ils(t, r_vec, ks, Cs)
        @test isa(g, Matrix)
        @test size(g) == (length(t), length(r_vec))
        for col in axes(g, 2)
            @test all(diff(g[:, col]) .> 0)   # increasing with time
        end
        for row in axes(g, 1)
            @test all(diff(g[row, :]) .< 0)   # decreasing with radius
        end
    end

    @testset "Vector overload agrees with scalar at each point" begin
        t = [3600.0, 7200.0, 14400.0]
        g_vec = ils(t, rb, ks, Cs)
        for i in eachindex(t)
            @test isapprox(g_vec[i], ils(t[i], rb, ks, Cs), rtol=1e-12)
        end
    end

end

# ICS — Infinite Cylindrical Source
@testset "ICS — Infinite Cylindrical Source" begin

    @testset "Scalar time and radius" begin
        g = ics(3600.0, rb, rb, ks, Cs)
        @test g > 0
        @test isfinite(g)
        @test isa(g, AbstractFloat)
    end

    @testset "Integer time is promoted to Float" begin
        @test isa(ics(3600, rb, rb, ks, Cs), AbstractFloat)
    end

    @testset "Vector time → vector, monotonically increasing" begin
        t = 3600.0 .* exp10.(range(0, log10(8760), length=20))
        g = ics(t, rb, rb, ks, Cs)
        @test isa(g, Vector)
        @test length(g) == length(t)
        @test all(isfinite, g)
        @test all(diff(g) .> 0)
    end

    @testset "Vector radius → vector, monotonically decreasing" begin
        r_vec = [0.076, 0.2, 0.5]
        g = ics(3600.0, r_vec, rb, ks, Cs)
        @test isa(g, Vector)
        @test all(isfinite, g)
        @test all(diff(g) .< 0)
    end

    @testset "Vector time × vector radius → matrix" begin
        t = 3600.0 .* exp10.(range(0, log10(8760), length=6))
        r_vec = [0.076, 0.2, 0.5]
        g = ics(t, r_vec, rb, ks, Cs)
        @test isa(g, Matrix)
        @test size(g) == (length(t), length(r_vec))
        for col in axes(g, 2)
            @test all(diff(g[:, col]) .> 0)
        end
        for row in axes(g, 1)
            @test all(diff(g[row, :]) .< 0)
        end
    end

    @testset "Vector overload agrees with scalar at each point" begin
        t = [3600.0, 7200.0, 14400.0]
        g_vec = ics(t, rb, rb, ks, Cs)
        for i in eachindex(t)
            @test isapprox(g_vec[i], ics(t[i], rb, rb, ks, Cs), rtol=1e-12)
        end
    end

    @testset "ICS ≈ ILS at the borehole wall for long times (within 5%)" begin
        t_long = 1e8
        @test isapprox(ics(t_long, rb, rb, ks, Cs), ils(t_long, rb, ks, Cs), rtol=0.05)
    end

    @testset "Different cylinder radii give valid g-functions" begin
        t = 3600.0
        for rc_test in [0.05, 0.076, 0.15]
            g = ics(t, 0.1, rc_test, ks, Cs)
            @test g > 0 && isfinite(g)
        end
    end
end

# FLS — Finite Line Source
@testset "FLS — Finite Line Source" begin

    @testset "Scalar time and radius" begin
        g = fls(3600.0, rb, H, D, ks, Cs)
        @test g > 0
        @test isfinite(g)
        @test isa(g, AbstractFloat)
    end

    @testset "Integer time is promoted to Float" begin
        @test isa(fls(3600, rb, H, D, ks, Cs), AbstractFloat)
    end

    @testset "Vector time → vector, monotonically increasing" begin
        t = 3600.0 .* exp10.(range(0, log10(8760 * 5), length=20))
        g = fls(t, rb, H, D, ks, Cs)
        @test isa(g, Vector)
        @test length(g) == length(t)
        @test all(isfinite, g)
        @test all(diff(g) .> 0)
    end

    @testset "Vector radius → vector, monotonically decreasing" begin
        r_vec = [0.076, 0.5, 2.0, 5.0]
        g = fls(3600.0 * 8760, r_vec, H, D, ks, Cs)
        @test isa(g, Vector)
        @test all(isfinite, g)
        @test all(diff(g) .< 0)
    end

    @testset "Vector time × vector radius → matrix" begin
        t = 3600.0 .* exp10.(range(0, 3, length=10))
        r_vec = [0.076, 1.0, 5.0]
        g = fls(t, r_vec, H, D, ks, Cs)
        @test isa(g, Matrix)
        @test size(g) == (length(t), length(r_vec))
    end

    @testset "Vector overload agrees with scalar at each point" begin
        t = [3600.0, 7200.0, 14400.0]
        g_vec = fls(t, rb, H, D, ks, Cs)
        for i in eachindex(t)
            @test isapprox(g_vec[i], fls(t[i], rb, H, D, ks, Cs), rtol=1e-12)
        end
    end

    @testset "FLS < ILS at long times (finite depth saturation)" begin
        # FLS saturates due to finite borehole depth; ILS keeps growing indefinitely
        t_long = 3600.0 * 8760 * 5
        @test fls(t_long, rb, H, D, ks, Cs) < ils(t_long, rb, ks, Cs)
    end

    @testset "Shallower borehole gives lower g-function" begin
        t_1yr = 3600.0 * 8760
        @test fls(t_1yr, rb, 75.0, D, ks, Cs) < fls(t_1yr, rb, H, D, ks, Cs)
    end

    @testset "Segment kernel reduces to classical FLS when H1=H2, D1=D2" begin
        t = 3600.0 .* exp10.(range(0, log10(8760 * 5), length=15))
        @test isapprox(fls(t, rb, H, D, H, D, ks, Cs), fls(t, rb, H, D, ks, Cs), rtol=1e-6)
    end

    @testset "Segment kernel reciprocity  H2·g(1→2) = H1·g(2→1)" begin
        t = 3600.0 * 8760
        H1, D1, H2, D2, r = 50.0, 4.0, 30.0, 80.0, 5.0
        g12 = fls(t, r, H1, D1, H2, D2, ks, Cs)
        g21 = fls(t, r, H2, D2, H1, D1, ks, Cs)
        @test isapprox(H2 * g12, H1 * g21, rtol=1e-6)
    end

    @testset "Segment kernel vector overload agrees with scalar" begin
        t = [3600.0, 7200.0, 14400.0]
        H1, D1, H2, D2 = 50.0, 4.0, 50.0, 54.0
        g_vec = fls(t, rb, H1, D1, H2, D2, ks, Cs)
        for i in eachindex(t)
            @test isapprox(g_vec[i], fls(t[i], rb, H1, D1, H2, D2, ks, Cs), rtol=1e-12)
        end
    end

end

# MILS — Moving Infinite Line Source
@testset "MILS — Moving Infinite Line Source" begin

    @testset "Scalar time at borehole wall (r = rb)" begin
        g = mils(3600.0, rb, 0.0, rb, ks, Cs, Cf, vD)
        @test g > 0
        @test isfinite(g)
        @test isa(g, AbstractFloat)
    end

    @testset "Scalar time at offset point (r = 5, θ = 0)" begin
        g = mils(3600.0*8760, 5.0, 0.0, rb, ks, Cs, Cf, vD)
        @test g > 0
        @test isfinite(g)
    end

    @testset "Integer time is promoted to Float" begin
        @test isa(mils(3600, rb, 0.0, rb, ks, Cs, Cf, vD), AbstractFloat)
    end

    @testset "Vector time → vector, monotonically increasing" begin
        t = 3600.0 .* exp10.(range(0, log10(8760), length=8))
        g = mils(t, rb, 0.0, rb, ks, Cs, Cf, vD)
        @test isa(g, Vector)
        @test length(g) == length(t)
        @test all(isfinite, g)
        @test all(diff(g) .> 0)
    end

    @testset "Vector overload agrees with scalar at each point" begin
        t = [3600.0, 7200.0, 14400.0]
        g_vec = mils(t, rb, 0.0, rb, ks, Cs, Cf, vD)
        for i in eachindex(t)
            @test isapprox(g_vec[i], mils(t[i], rb, 0.0, rb, ks, Cs, Cf, vD), rtol=1e-12)
        end
    end

    @testset "Downstream (θ = 0) warmer than upstream (θ = 180) at equal distance" begin
        # Groundwater flows in +x; plume extends downstream
        t_1yr = 3600.0 * 8760
        r_off = 5.0
        g_down = mils(t_1yr, r_off,   0.0, rb, ks, Cs, Cf, vD)
        g_up   = mils(t_1yr, r_off, 180.0, rb, ks, Cs, Cf, vD)
        @test g_down > g_up
    end

    @testset "g decreases along downstream x-axis with distance" begin
        t_1yr = 3600.0 * 8760
        g1 = mils(t_1yr,  2.0, 0.0, rb, ks, Cs, Cf, vD)
        g2 = mils(t_1yr,  5.0, 0.0, rb, ks, Cs, Cf, vD)
        g3 = mils(t_1yr, 10.0, 0.0, rb, ks, Cs, Cf, vD)
        @test g1 > g2 > g3
    end

    @testset "Very low vD → approaches ILS (within 5%)" begin
        t_1yr  = 3600.0 * 8760
        vD_min = 1e-12
        g_mils = mils(t_1yr, rb, 0.0, rb, ks, Cs, Cf, vD_min)
        g_ils  = ils(t_1yr, rb, ks, Cs)
        @test isapprox(g_mils, g_ils, rtol=0.05)
    end

    @testset "Self response is the circumferential mean at the wall (r = rb)" begin
        # The self term must be the mean borehole-wall temperature, i.e. the θ-average of the
        # directional field at r = rb — NOT the borehole-centre value. Guard at high vD, where
        # the two differ by the factor I₀(rb·vT/2α) ≫ 1 (regression against the old r = 0 self).
        t_1yr = 3600.0 * 8760
        vD_hi = 1e-4
        self  = mils(t_1yr, rb, 0.0, rb, ks, Cs, Cf, vD_hi)
        nθ    = 2000
        rw    = rb * (1 + 1e-4)                        # just outside the wall (directional branch)
        mean  = sum(mils(t_1yr, rw, θ, rb, ks, Cs, Cf, vD_hi)
                    for θ in range(0, 360, length = nθ + 1)[1:end-1]) / nθ
        @test isapprox(self, mean, rtol = 1e-3)
    end

    @testset "Borefield matrix overload (nb×nb) — shape and finiteness" begin
        xy22 = borefield(:rectangle, 2, 2, B)
        nb   = size(xy22, 1)
        r, θ = borefield_geometry(xy22, rb)
        g2D  = mils(3600.0, r, θ, rb, ks, Cs, Cf, vD)
        @test isa(g2D, Matrix)
        @test size(g2D) == (nb, nb)
        @test all(isfinite, g2D)
    end

    @testset "3D borefield overload (nt×nb×nb) — shape and finiteness" begin
        xy22 = borefield(:rectangle, 2, 2, B)
        nb   = size(xy22, 1)
        t    = [3600.0, 7200.0, 14400.0]
        r, θ = borefield_geometry(xy22, rb)
        g3D  = mils(t, r, θ, rb, ks, Cs, Cf, vD)
        @test ndims(g3D) == 3
        @test size(g3D) == (length(t), nb, nb)
        @test all(isfinite, g3D)
    end

end

# MFLS — Moving Finite Line Source
@testset "MFLS — Moving Finite Line Source" begin

    @testset "Scalar time at borehole wall (r = rb)" begin
        g = mfls(3600.0, rb, 0.0, H, rb, D, ks, Cs, Cf, vD)
        @test g > 0
        @test isfinite(g)
        @test isa(g, AbstractFloat)
    end

    @testset "Scalar time at offset point (r = 5, θ = 0)" begin
        g = mfls(3600.0*8760, 5.0, 0.0, H, rb, D, ks, Cs, Cf, vD)
        @test g > 0
        @test isfinite(g)
    end

    @testset "Integer time is promoted to Float" begin
        @test isa(mfls(3600, rb, 0.0, H, rb, D, ks, Cs, Cf, vD), AbstractFloat)
    end

    @testset "Vector time → vector, monotonically increasing" begin
        t = 3600.0 .* exp10.(range(0, 3, length=8))
        g = mfls(t, rb, 0.0, H, rb, D, ks, Cs, Cf, vD)
        @test isa(g, Vector)
        @test length(g) == length(t)
        @test all(isfinite, g)
        @test all(diff(g) .> 0)
    end

    @testset "Vector overload agrees with scalar at each point" begin
        t = [3600.0, 7200.0]
        g_vec = mfls(t, rb, 0.0, H, rb, D, ks, Cs, Cf, vD)
        for i in eachindex(t)
            @test isapprox(g_vec[i], mfls(t[i], rb, 0.0, H, rb, D, ks, Cs, Cf, vD), rtol=1e-12)
        end
    end

    @testset "Downstream (θ = 0) warmer than upstream (θ = 180) at equal distance" begin
        t_1yr = 3600.0 * 8760
        r_off = 5.0
        g_down = mfls(t_1yr, r_off,   0.0, H, rb, D, ks, Cs, Cf, vD)
        g_up   = mfls(t_1yr, r_off, 180.0, H, rb, D, ks, Cs, Cf, vD)
        @test g_down > g_up
    end

    @testset "MFLS < MILS at long times (finite depth effect)" begin
        # MFLS saturates due to finite borehole depth; MILS is infinite
        t_long = 3600.0 * 8760 * 5
        g_mfls = mfls(t_long, rb, 0.0, H, rb, D, ks, Cs, Cf, vD)
        g_mils = mils(t_long, rb, 0.0, rb, ks, Cs, Cf, vD)
        @test g_mfls < g_mils
    end

    @testset "Very low vD → approaches FLS (within 5%)" begin
        t_1yr  = 3600.0 * 8760
        vD_min = 1e-12
        g_mfls = mfls(t_1yr, rb, 0.0, H, rb, D, ks, Cs, Cf, vD_min)
        g_fls  = fls(t_1yr, rb, H, D, ks, Cs)
        @test isapprox(g_mfls, g_fls, rtol=0.05)
    end

    @testset "Self response is the circumferential mean at the wall (r = rb)" begin
        # As for MILS: the self term must be the mean borehole-wall temperature, i.e. the θ-average
        # of the directional field at r = rb — NOT the borehole-centre value. Guard at high vD,
        # where the two differ by the factor I₀(rb·U/2α) ≫ 1 (regression against the old r = 0 self).
        t_1yr = 3600.0 * 8760
        vD_hi = 1e-4
        self  = mfls(t_1yr, rb, 0.0, H, rb, D, ks, Cs, Cf, vD_hi)
        nθ    = 2000
        rw    = rb * (1 + 1e-4)                        # just outside the wall (directional branch)
        mean  = sum(mfls(t_1yr, rw, θ, H, rb, D, ks, Cs, Cf, vD_hi)
                    for θ in range(0, 360, length = nθ + 1)[1:end-1]) / nθ
        @test isapprox(self, mean, rtol = 1e-3)
    end
end

# Spatial superposition
@testset "borefield_geometry" begin
    xy = borefield(:rectangle, 3, 2, 5.0)
    r, θ = borefield_geometry(xy, rb)
    nb = size(xy, 1)
    @test size(r) == (nb, nb) && size(θ) == (nb, nb)
    @test all(diag(r) .== rb)                      # self-distance on the diagonal equals rb
    @test all(diag(θ) .== 0)                       # angle 0 on the diagonal
    # downstream (+x) receiver has θ = 0°, upstream θ = 180°
    xline = borefield(:line, 3, 5.0)               # (0,0),(5,0),(10,0)
    _, θl, = borefield_geometry(xline, rb)
    @test θl[3, 1] ≈ 0.0   atol=1e-9               # receiver 3 downstream of source 1
    @test θl[1, 3] ≈ 180.0 atol=1e-9               # receiver 1 upstream of source 3
end

# Borefield layout functions
@testset "Borefield layout functions" begin

    @testset "borefield_rectangle — borehole count and spacing" begin
        nx, ny = 3, 4
        xy = borefield(:rectangle, nx, ny, B)
        @test size(xy) == (nx * ny, 2)
        @test length(unique(round.(xy[:, 1], digits=8))) == nx
        @test length(unique(round.(xy[:, 2], digits=8))) == ny
    end

    @testset "borefield_rectangle — anisotropic spacing" begin
        Bx, By = 4.0, 7.0
        xy = borefield(:rectangle, 3, 4, Bx, By)
        @test size(xy) == (12, 2)
        xs = sort(unique(round.(xy[:, 1], digits=8)))
        ys = sort(unique(round.(xy[:, 2], digits=8)))
        @test isapprox(xs[2] - xs[1], Bx, rtol=1e-10)
        @test isapprox(ys[2] - ys[1], By, rtol=1e-10)
    end

    @testset "borefield_line — count and y = 0" begin
        n = 6
        xy = borefield(:line, n, B)
        @test size(xy) == (n, 2)
        @test all(xy[:, 2] .== 0.0)
        xs = sort(xy[:, 1])
        @test isapprox(xs[2] - xs[1], B, rtol=1e-10)
    end

    @testset "borefield_circle — count and all points at radius R" begin
        nb = 12
        R  = 20.0
        xy = borefield(:circle, nb, R)
        @test size(xy) == (nb, 2)
        r_vals = sqrt.(xy[:, 1].^2 .+ xy[:, 2].^2)
        @test all(isapprox.(r_vals, R, rtol=1e-10))
    end

    @testset "borefield_L — borehole count" begin
        n1, n2 = 6, 5
        xy = borefield(:L, n1, n2, B)
        @test size(xy, 1) == n1 + n2 - 1
        @test size(xy, 2) == 2
    end

    @testset "borefield_L — anisotropic spacing" begin
        n1, n2 = 6, 5
        xy = borefield(:L, n1, n2, 4.0, 6.0)
        @test size(xy, 1) == n1 + n2 - 1
    end

    @testset "borefield_U — borehole count" begin
        nx, ny = 5, 4
        xy = borefield(:U, nx, ny, B)
        @test size(xy, 1) == nx + 2 * (ny - 1)
        @test size(xy, 2) == 2
    end

    @testset "borefield_open_rectangle — borehole count" begin
        nx, ny = 6, 4
        xy = borefield(:open_rectangle, nx, ny, B)
        @test size(xy, 1) == 2 * (nx + ny - 2)
        @test size(xy, 2) == 2
    end

    @testset "borefield wrapper — unknown shape throws ArgumentError" begin
        @test_throws ArgumentError borefield(:hexagon, 5, B)
    end

    @testset "borefield_geometry — output shape and diagonal (rb, 0°)" begin
        xy = borefield(:rectangle, 3, 3, B)
        nb = size(xy, 1)
        r, θ = borefield_geometry(xy, rb)
        @test size(r) == (nb, nb)
        @test size(θ) == (nb, nb)
        for i in 1:nb
            @test r[i, i] == rb
            @test θ[i, i] == 0
        end
    end

    @testset "borefield_geometry — off-diagonal distances exceed rb" begin
        xy = borefield(:rectangle, 2, 2, B)
        r, = borefield_geometry(xy, rb)
        nb = size(r, 1)
        for i in 1:nb, j in 1:nb
            i != j && @test r[i, j] > rb
        end
    end
end

@testset "Borefield rotation (angle keyword)" begin
    # Counter-clockwise rotation about the origin, in degrees
    @test borefield_line(3, B; angle = 90) ≈ [0.0 0.0; 0.0 B; 0.0 2B] atol = 1e-12
    @test borefield(:rectangle, 2, 2, B; angle = 30) ≈ borefield_rectangle(2, 2, B; angle = 30)
    # Distances unchanged, angles to the flow changed
    r0, θ0 = borefield_geometry(borefield_rectangle(3, 2, B), rb)
    r1, θ1 = borefield_geometry(borefield_rectangle(3, 2, B; angle = 30), rb)
    @test r1 ≈ r0
    @test !(θ1 ≈ θ0)
end

# Response arrays — gmatrix
@testset "gmatrix — response arrays" begin

    t    = 2_592_000.0 .* (1:12)
    xy1  = [0.0 0.0]
    xy22 = borefield(:rectangle, 2, 2, B)
    m    = FLSModel(H, D, ks, Cs)

    @testset "Single borehole — the kernel at the wall" begin
        @test gmatrix(t, rb, xy1, ILSModel(ks, Cs)) ≈ ils(t, rb, ks, Cs)
        @test gmatrix(t, rb, xy1, ICSModel(ks, Cs)) ≈ ics(t, rb, rb, ks, Cs)
        @test gmatrix(t, rb, xy1, m) ≈ fls(t, rb, H, D, ks, Cs)
        @test gmatrix(t, rb, xy1, MILSModel(rb, ks, Cs, Cf, vD)) ≈ mils(t, rb, 0.0, rb, ks, Cs, Cf, vD)
        @test gmatrix(t, rb, xy1, MFLSModel(H, rb, D, ks, Cs, Cf, vD)) ≈
            mfls(t, rb, 0.0, H, rb, D, ks, Cs, Cf, vD)
    end

    @testset "Borefield — self response on the diagonal, symmetric for an isotropic model" begin
        g = gmatrix(t, rb, xy22, m)
        @test size(g) == (length(t), 4, 4)
        @test all(g[:, i, i] ≈ fls(t, rb, H, D, ks, Cs) for i in 1:4)
        @test all(g[:, i, j] ≈ g[:, j, i] for i in 1:4, j in 1:4)
    end

    @testset "Moving model needs its own radius" begin
        @test_throws ArgumentError gmatrix(t, 2rb, xy22, MILSModel(rb, ks, Cs, Cf, vD))
    end

    @testset "Segments — one segment is the borehole array" begin
        @test gmatrix(t, rb, xy22, m, 1) ≈ gmatrix(t, rb, xy22, m) rtol = 1e-6
        @test_throws ArgumentError gmatrix(t, rb, xy22, m, 0)
    end

    @testset "Segments — array built by hand" begin
        nseg = 6; nb = size(xy22, 1); Hs = H / nseg
        Dseg = [D + (k - 1) * Hs for k in 1:nseg]
        r, = borefield_geometry(xy22, rb); rmat = max.(r, rb)
        seg = [(bh, k) for bh in 1:nb for k in 1:nseg]; nS = nb * nseg
        g3D = zeros(length(t), nS, nS)
        for q in 1:nS, p in 1:nS
            ib, is = seg[p]; jb, js = seg[q]
            g3D[:, p, q] = fls(t, rmat[ib, jb], Hs, Dseg[js], Hs, Dseg[is], ks, Cs)
        end
        @test gmatrix(t, rb, xy22, m, nseg) ≈ g3D rtol = 1e-12
    end

    @testset "Directional (r,θ) dedup build equals naive per-pair build" begin
        local m = MFLSModel(H, rb, D, ks, Cs, Cf, vD)
        xy33 = borefield(:rectangle, 3, 3, B)
        r, θ = borefield_geometry(xy33, m.rb)
        g_dedup = mfls(t, r, θ, m.H, m.rb, m.D, m.ks, m.Cs, m.Cf, m.vD)
        nb = size(xy33, 1)
        g_naive = [mfls(t[k], r[i,j], θ[i,j], m.H, m.rb, m.D, m.ks, m.Cs, m.Cf, m.vD)
                   for k in eachindex(t), i in 1:nb, j in 1:nb]
        @test g_dedup == g_naive
    end

    @testset "Rotation — isotropic array unchanged, moving array changed" begin
        xy = borefield_rectangle(3, 2, B)
        @test gmatrix(t, rb, borefield_rectangle(3, 2, B; angle = 30), m) ≈ gmatrix(t, rb, xy, m)
        mm = MILSModel(rb, ks, Cs, Cf, vD)
        @test !(gmatrix(t, rb, borefield_rectangle(3, 2, B; angle = 90), mm) ≈ gmatrix(t, rb, xy, mm))
    end
end
