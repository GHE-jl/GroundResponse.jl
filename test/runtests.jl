using Test
using GroundResponse

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

end

# MILS — Moving Infinite Line Source
@testset "MILS — Moving Infinite Line Source" begin

    @testset "Scalar time at borehole wall [0, 0]" begin
        g = mils(3600.0, [0.0, 0.0], rb, ks, Cs, Cf, vD)
        @test g > 0
        @test isfinite(g)
        @test isa(g, AbstractFloat)
    end

    @testset "Scalar time at offset point [5, 0]" begin
        g = mils(3600.0*8760, [5.0, 0.0], rb, ks, Cs, Cf, vD)
        @test g > 0
        @test isfinite(g)
    end

    @testset "Integer time is promoted to Float" begin
        @test isa(mils(3600, [0.0, 0.0], rb, ks, Cs, Cf, vD), AbstractFloat)
    end

    @testset "Vector time → vector, monotonically increasing" begin
        t = 3600.0 .* exp10.(range(0, log10(8760), length=8))
        g = mils(t, [0.0, 0.0], rb, ks, Cs, Cf, vD)
        @test isa(g, Vector)
        @test length(g) == length(t)
        @test all(isfinite, g)
        @test all(diff(g) .> 0)
    end

    @testset "Vector overload agrees with scalar at each point" begin
        t = [3600.0, 7200.0, 14400.0]
        g_vec = mils(t, [0.0, 0.0], rb, ks, Cs, Cf, vD)
        for i in eachindex(t)
            @test isapprox(g_vec[i], mils(t[i], [0.0, 0.0], rb, ks, Cs, Cf, vD), rtol=1e-12)
        end
    end

    @testset "Downstream (+x) warmer than upstream (−x) at equal distance" begin
        # Groundwater flows in +x; plume extends downstream
        t_1yr = 3600.0 * 8760
        r_off = 5.0
        g_down = mils(t_1yr, [ r_off, 0.0], rb, ks, Cs, Cf, vD)
        g_up   = mils(t_1yr, [-r_off, 0.0], rb, ks, Cs, Cf, vD)
        @test g_down > g_up
    end

    @testset "g decreases along downstream x-axis with distance" begin
        t_1yr = 3600.0 * 8760
        g1 = mils(t_1yr, [ 2.0, 0.0], rb, ks, Cs, Cf, vD)
        g2 = mils(t_1yr, [ 5.0, 0.0], rb, ks, Cs, Cf, vD)
        g3 = mils(t_1yr, [10.0, 0.0], rb, ks, Cs, Cf, vD)
        @test g1 > g2 > g3
    end

    @testset "Very low vD → approaches ILS (within 5%)" begin
        t_1yr  = 3600.0 * 8760
        vD_min = 1e-12
        g_mils = mils(t_1yr, [0.0, 0.0], rb, ks, Cs, Cf, vD_min)
        g_ils  = ils(t_1yr, rb, ks, Cs)
        @test isapprox(g_mils, g_ils, rtol=0.05)
    end

    @testset "Borefield matrix overload (nb×nb) — shape and finiteness" begin
        xy22 = borefield(:rectangle, 2, 2, B)
        nb   = size(xy22, 1)
        g2D  = mils(3600.0, xy22, rb, ks, Cs, Cf, vD)
        @test isa(g2D, Matrix)
        @test size(g2D) == (nb, nb)
        @test all(isfinite, g2D)
    end

    @testset "3D borefield overload (nt×nb×nb) — shape and finiteness" begin
        xy22 = borefield(:rectangle, 2, 2, B)
        nb   = size(xy22, 1)
        t    = [3600.0, 7200.0, 14400.0]
        g3D  = mils(t, xy22, rb, ks, Cs, Cf, vD)
        @test ndims(g3D) == 3
        @test size(g3D) == (length(t), nb, nb)
        @test all(isfinite, g3D)
    end

end

# MFLS — Moving Finite Line Source
@testset "MFLS — Moving Finite Line Source" begin

    @testset "Scalar time at borehole wall [0, 0]" begin
        g = mfls(3600.0, [0.0, 0.0], H, rb, D, ks, Cs, Cf, vD)
        @test g > 0
        @test isfinite(g)
        @test isa(g, AbstractFloat)
    end

    @testset "Scalar time at offset point [5, 0]" begin
        g = mfls(3600.0*8760, [5.0, 0.0], H, rb, D, ks, Cs, Cf, vD)
        @test g > 0
        @test isfinite(g)
    end

    @testset "Integer time is promoted to Float" begin
        @test isa(mfls(3600, [0.0, 0.0], H, rb, D, ks, Cs, Cf, vD), AbstractFloat)
    end

    @testset "Vector time → vector, monotonically increasing" begin
        t = 3600.0 .* exp10.(range(0, 3, length=8))
        g = mfls(t, [0.0, 0.0], H, rb, D, ks, Cs, Cf, vD)
        @test isa(g, Vector)
        @test length(g) == length(t)
        @test all(isfinite, g)
        @test all(diff(g) .> 0)
    end

    @testset "Vector overload agrees with scalar at each point" begin
        t = [3600.0, 7200.0]
        g_vec = mfls(t, [0.0, 0.0], H, rb, D, ks, Cs, Cf, vD)
        for i in eachindex(t)
            @test isapprox(g_vec[i], mfls(t[i], [0.0, 0.0], H, rb, D, ks, Cs, Cf, vD), rtol=1e-12)
        end
    end

    @testset "Downstream (+x) warmer than upstream (−x) at equal distance" begin
        t_1yr = 3600.0 * 8760
        r_off = 5.0
        g_down = mfls(t_1yr, [ r_off, 0.0], H, rb, D, ks, Cs, Cf, vD)
        g_up   = mfls(t_1yr, [-r_off, 0.0], H, rb, D, ks, Cs, Cf, vD)
        @test g_down > g_up
    end

    @testset "MFLS < MILS at long times (finite depth effect)" begin
        # MFLS saturates due to finite borehole depth; MILS is infinite
        t_long = 3600.0 * 8760 * 5
        g_mfls = mfls(t_long, [0.0, 0.0], H, rb, D, ks, Cs, Cf, vD)
        g_mils = mils(t_long, [0.0, 0.0], rb, ks, Cs, Cf, vD)
        @test g_mfls < g_mils
    end

    @testset "Very low vD → approaches FLS (within 5%)" begin
        t_1yr  = 3600.0 * 8760
        vD_min = 1e-12
        g_mfls = mfls(t_1yr, [0.0, 0.0], H, rb, D, ks, Cs, Cf, vD_min)
        g_fls  = fls(t_1yr, rb, H, D, ks, Cs)
        @test isapprox(g_mfls, g_fls, rtol=0.05)
    end

end

# Spatial superposition
@testset "Spatial superposition" begin

    m_fls = FLSModel(H, D, ks, Cs)
    t     = 3600.0 .* exp10.(range(0, log10(8760 * 10), length=20))
    xy22  = borefield(:rectangle, 2, 2, B)
    xy33  = borefield(:rectangle, 3, 3, B)

    @testset "successive_flux and bloc_matrix agree — 2×2 FLS" begin
        g_sf = successive_flux(t, rb, xy22, m_fls)
        g_bm = bloc_matrix(t, rb, xy22, m_fls)
        @test isapprox(g_sf, g_bm, rtol=1e-4)
    end

    @testset "successive_flux and bloc_matrix agree — 3×3 FLS" begin
        g_sf = successive_flux(t, rb, xy33, m_fls)
        g_bm = bloc_matrix(t, rb, xy33, m_fls)
        @test isapprox(g_sf, g_bm, rtol=1e-2)
    end

    @testset "Borefield g-function smaller than single-borehole response since more efficient" begin
        g_1   = ils(t, rb, ks, Cs)
        g_4   = successive_flux(t, rb, xy22, ILSModel(ks, Cs))
        @test all(g_4 .< g_1)
    end

    @testset "Larger borefield gives smaller g-function" begin
        g22 = successive_flux(t, rb, xy22, m_fls)
        g33 = successive_flux(t, rb, xy33, m_fls)
        @test all(g33 .< g22)
    end

    @testset "successive_flux output is monotonically increasing" begin
        g = successive_flux(t, rb, xy22, m_fls)
        @test all(diff(g) .> 0)
    end

    @testset "Low-level 3D array overload matches high-level — successive_flux" begin
        r, = borefield_radius(xy22, rb)
        g3D  = fls(t, r, H, D, ks, Cs)
        g_ll = successive_flux(g3D)
        g_hl = successive_flux(t, rb, xy22, m_fls)
        @test isapprox(g_ll, g_hl, rtol=1e-10)
    end

    @testset "Low-level 3D array overload matches high-level — bloc_matrix" begin
        r, = borefield_radius(xy22, rb)
        g3D  = fls(t, r, H, D, ks, Cs)
        g_ll = bloc_matrix(g3D)
        g_hl = bloc_matrix(t, rb, xy22, m_fls)
        @test isapprox(g_ll, g_hl, rtol=1e-10)
    end

    @testset "successive_flux works with MILSModel" begin
        m_mils = MILSModel(rb, ks, Cs, Cf, vD)
        g = successive_flux(t, rb, xy22, m_mils)
        @test length(g) == length(t)
        @test all(isfinite, g)
        @test all(g .> 0)
    end
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

    @testset "borefield_radius — output shape and diagonal equals rb" begin
        xy = borefield(:rectangle, 3, 3, B)
        nb = size(xy, 1)
        r, rᵥ, rᵤ, rᵢ, θ, nb_out = borefield_radius(xy, rb)
        @test nb_out == nb
        @test size(r) == (nb, nb)
        @test length(rᵥ) == nb^2
        @test length(θ)  == nb
        for i in 1:nb
            @test r[i, i] == rb
        end
    end

    @testset "borefield_radius — off-diagonal entries exceed rb" begin
        xy = borefield(:rectangle, 2, 2, B)
        r, = borefield_radius(xy, rb)
        nb = size(r, 1)
        for i in 1:nb, j in 1:nb
            i != j && @test r[i, j] > rb
        end
    end
end

# High-level interface — ground_response
@testset "ground_response — high-level interface" begin

    t    = 3600.0 .* exp10.(range(0, log10(8760 * 5), length=30))
    xy1  = [0.0 0.0]   # single borehole (1×2 matrix)
    xy22 = borefield(:rectangle, 2, 2, B)

    @testset "Single borehole — ILSModel matches ils directly" begin
        m = ILSModel(ks, Cs)
        @test isapprox(ground_response(t, rb, xy1, m), ils(t, rb, ks, Cs), rtol=1e-12)
    end

    @testset "Single borehole — ICSModel matches ics directly" begin
        m = ICSModel(rb, ks, Cs)
        @test isapprox(ground_response(t, rb, xy1, m), ics(t, rb, rb, ks, Cs), rtol=1e-12)
    end

    @testset "Single borehole — FLSModel matches fls directly" begin
        m = FLSModel(H, D, ks, Cs)
        @test isapprox(ground_response(t, rb, xy1, m), fls(t, rb, H, D, ks, Cs), rtol=1e-12)
    end

    @testset "Single borehole — MILSModel matches mils at [0, 0]" begin
        m = MILSModel(rb, ks, Cs, Cf, vD)
        @test isapprox(
            ground_response(t, rb, xy1, m),
            mils(t, [0.0, 0.0], rb, ks, Cs, Cf, vD),
            rtol=1e-12)
    end

    @testset "Single borehole — MFLSModel matches mfls at [0, 0]" begin
        m = MFLSModel(H, rb, D, ks, Cs, Cf, vD)
        @test isapprox(
            ground_response(t, rb, xy1, m),
            mfls(t, [0.0, 0.0], H, rb, D, ks, Cs, Cf, vD),
            rtol=1e-12)
    end

    @testset "Borefield dispatch matches successive_flux directly" begin
        m = FLSModel(H, D, ks, Cs)
        @test isapprox(
            ground_response(t, rb, xy22, m),
            successive_flux(t, rb, xy22, m),
            rtol=1e-12)
    end

    @testset "Borefield output is monotonically increasing" begin
        m = FLSModel(H, D, ks, Cs)
        g = ground_response(t, rb, xy22, m)
        @test all(diff(g) .> 0)
    end

    @testset "Borefield response smaller than single-borehole response" begin
        m = FLSModel(H, D, ks, Cs)
        g1  = ground_response(t, rb, xy1,  m)
        g22 = ground_response(t, rb, xy22, m)
        @test all(g22 .< g1)
    end
end
