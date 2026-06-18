"""
    test_infinite_line_source.jl

Tests for the Infinite Line Source (ILS) model.
"""

@testset "Infinite Line Source (ILS)" begin
    # Test parameters
    r = 0.076               # Borehole radius [m]
    ks = 3.0                # Ground thermal conductivity [W/mK]
    Cs = 2e6                # Ground volumetric specific heat [J/m³K]
    
    @testset "Single time, single radius" begin
        t = 3600.0
        g = ils(t, r, ks, Cs)
        
        @test g > 0.0
        @test isfinite(g)
        @test g < 1.0
    end
    
    @testset "Vector time, single radius" begin
        t = 60:60:3600
        g = ils(t, r, ks, Cs)
        
        @test isa(g, Vector)
        @test length(g) == length(t)
        @test all(g .> 0.0)
        @test all(isfinite.(g))
        @test all(diff(g) .> 0.0)  # Monotonically increasing
    end
    
    @testset "Single time, vector radius" begin
        t = 3600.0
        r_vec = [0.05, 0.076, 0.1, 0.15]
        g = ils(t, r_vec, ks, Cs)
        
        @test isa(g, Vector)
        @test length(g) == length(r_vec)
        @test all(g .> 0.0)
        @test all(isfinite.(g))
        @test all(diff(g) .< 0.0)  # Decreasing with radius
    end
    
    @testset "Vector time, vector radius (2D matrix)" begin
        t = 60:60:3600
        r_vec = [0.05, 0.076, 0.1]
        g = ils(t, r_vec, ks, Cs)
        
        @test isa(g, Matrix)
        @test size(g) == (length(t), length(r_vec))
        @test all(g .> 0.0)
        @test all(isfinite.(g))
        
        # Each column should increase with time
        for col in 1:size(g, 2)
            @test all(diff(g[:, col]) .> 0.0)
        end
        
        # Each row should decrease with radius
        for row in 1:size(g, 1)
            @test all(diff(g[row, :]) .< 0.0)
        end
    end
    
    @testset "Comparison with FLS" begin
        # ILS should be different from FLS (FLS includes buried depth)
        H = 150.0
        D = 4.0
        t = 3600.0
        
        g_ils = ils(t, r, ks, Cs)
        g_fls = fls(t, H, r, D, ks, Cs)
        
        @test g_ils ≠ g_fls
    end
    
    @testset "Type promotion" begin
        t_int = 3600
        g1 = ils(t_int, r, ks, Cs)
        @test isa(g1, AbstractFloat)
    end
end
