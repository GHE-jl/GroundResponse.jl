"""
    test_moving_finite_line_source.jl

Tests for the Moving Finite Line Source (MFLS) model.
"""

@testset "Moving Finite Line Source (MFLS)" begin
    # Test parameters
    H = 150.0               # Borehole depth [m]
    rb = 0.076              # Borehole radius [m]
    D = 4.0                 # Buried depth [m]
    ks = 3.0                # Ground thermal conductivity [W/mK]
    Cs = 2e6                # Ground volumetric specific heat [J/m³K]
    Cf = 4.2e6              # Groundwater volumetric specific heat [J/m³K]
    vD = 1e-6               # Darcy velocity [m/s]
    
    @testset "Single time, single coordinate" begin
        t = 3600.0
        xy = [0.0; 0.0]    # At borehole wall
        g = mfls(t, H, rb, D, xy, ks, Cs, Cf, vD)
        
        @test g > 0.0
        @test isfinite(g)
    end
    
    @testset "Vector time, single coordinate" begin
        t = 60:60:3600
        xy = [0.0; 0.0]
        g = mfls(t, H, rb, D, xy, ks, Cs, Cf, vD)
        
        @test isa(g, Vector)
        @test length(g) == length(t)
        @test all(g .> 0.0)
        @test all(isfinite.(g))
        @test all(diff(g) .> 0.0)  # Monotonically increasing
    end
    
    @testset "Single time, multiple coordinates" begin
        t = 3600.0
        xy = [0.0 5.0 10.0; 0.0 0.0 0.0]  # 3 different x positions
        g = mfls(t, H, rb, D, xy, ks, Cs, Cf, vD)
        
        @test isa(g, Vector)
        @test length(g) == size(xy, 2)
        @test all(g .> 0.0)
        @test all(isfinite.(g))
        
        # g-function should decrease with distance from borehole
        @test all(diff(g) .< 0.0)
    end
    
    @testset "Vector time, multiple coordinates" begin
        t = 60:60:3600
        xy = [0.0 5.0 10.0; 0.0 0.0 0.0]
        g = mfls(t, H, rb, D, xy, ks, Cs, Cf, vD)
        
        @test isa(g, Matrix)
        @test size(g) == (length(t), size(xy, 2))
        @test all(g .> 0.0)
        @test all(isfinite.(g))
        
        # Each column should increase with time
        for col in 1:size(g, 2)
            @test all(diff(g[:, col]) .> 0.0)
        end
        
        # Each row should decrease with distance
        for row in 1:size(g, 1)
            @test all(diff(g[row, :]) .< 0.0)
        end
    end
    
    @testset "Inside vs outside borehole" begin
        t = 3600.0
        xy_inside = [0.0; 0.0]     # Inside borehole
        xy_outside = [5.0; 0.0]    # Outside borehole
        
        g_inside = mfls(t, H, rb, D, xy_inside, ks, Cs, Cf, vD)
        g_outside = mfls(t, H, rb, D, xy_outside, ks, Cs, Cf, vD)
        
        @test g_inside > 0.0
        @test g_outside > 0.0
        @test isfinite(g_inside)
        @test isfinite(g_outside)
        
        # Inside should give larger g-function
        @test g_inside > g_outside
    end
    
    @testset "Different groundwater velocities" begin
        t = 3600.0
        xy = [0.0; 0.0]
        
        vD1 = 1e-12  # Very low velocity (nearly impervious)
        vD2 = 1e-6   # Low velocity
        vD3 = 1e-5   # Higher velocity
        
        g1 = mfls(t, H, rb, D, xy, ks, Cs, Cf, vD1)
        g2 = mfls(t, H, rb, D, xy, ks, Cs, Cf, vD2)
        g3 = mfls(t, H, rb, D, xy, ks, Cs, Cf, vD3)
        
        @test all([g1, g2, g3] .> 0.0)
        @test all(isfinite.([g1, g2, g3]))
    end
    
    @testset "Type promotion" begin
        t_int = 3600
        xy = [0.0; 0.0]
        g1 = mfls(t_int, H, rb, D, xy, ks, Cs, Cf, vD)
        @test isa(g1, AbstractFloat)
    end
end
