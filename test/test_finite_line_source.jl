"""
    test_finite_line_source.jl

Tests for the Finite Line Source (FLS) model.
"""

@testset "Finite Line Source (FLS)" begin
    # Test parameters - typical borehole configuration
    H = 150.0               # Borehole depth [m]
    r = 0.076               # Borehole radius [m]
    D = 4.0                 # Buried depth [m]
    ks = 3.0                # Ground thermal conductivity [W/mK]
    Cs = 2e6                # Ground volumetric specific heat [J/m³K]
    
    @testset "Single time, single radius" begin
        t = 3600.0          # 1 hour [s]
        g = fls(t, H, r, D, ks, Cs)
        
        # g-function should be positive and finite
        @test g > 0.0
        @test isfinite(g)
        # g-function should be small for short times (logarithmic growth)
        @test g < 1.0
    end
    
    @testset "Vector time, single radius" begin
        t = 60:60:3600      # Time vector from 60 to 3600 seconds
        g = fls(t, H, r, D, ks, Cs)
        
        # Output should be a vector with same length as input
        @test isa(g, Vector)
        @test length(g) == length(t)
        
        # All values should be positive and finite
        @test all(g .> 0.0)
        @test all(isfinite.(g))
        
        # g-function should be monotonically increasing with time
        @test all(diff(g) .> 0.0)
    end
    
    @testset "Single time, vector radius" begin
        t = 3600.0
        r_vec = [0.05, 0.076, 0.1, 0.15]
        g = fls(t, H, r_vec, D, ks, Cs)
        
        # Output should be a vector with same length as radius vector
        @test isa(g, Vector)
        @test length(g) == length(r_vec)
        
        # All values should be positive and finite
        @test all(g .> 0.0)
        @test all(isfinite.(g))
        
        # g-function should decrease with increasing radius (distance from borehole)
        @test all(diff(g) .< 0.0)
    end
    
    @testset "Vector time, vector radius (2D matrix)" begin
        t = 60:60:3600
        r_vec = [0.05, 0.076, 0.1]
        g = fls(t, H, r_vec, D, ks, Cs)
        
        # Output should be a matrix: rows = time steps, columns = radius values
        @test isa(g, Matrix)
        @test size(g) == (length(t), length(r_vec))
        
        # All values should be positive and finite
        @test all(g .> 0.0)
        @test all(isfinite.(g))
        
        # Each column should increase with time (fixed radius)
        for col in 1:size(g, 2)
            @test all(diff(g[:, col]) .> 0.0)
        end
        
        # Each row should decrease with radius (fixed time)
        for row in 1:size(g, 1)
            @test all(diff(g[row, :]) .< 0.0)
        end
    end
    
    @testset "Type promotion" begin
        # Test with different floating point types
        t_int = 3600
        g1 = fls(t_int, H, r, D, ks, Cs)
        @test isa(g1, AbstractFloat)
        
        # Test with mixed types
        t_float32 = Float32(3600.0)
        g2 = fls(t_float32, H, r, D, ks, Cs)
        @test isa(g2, AbstractFloat)
    end
    
    @testset "Edge cases" begin
        # Very short time (should give small g-value)
        t_short = 1.0
        g_short = fls(t_short, H, r, D, ks, Cs)
        @test g_short > 0.0
        @test isfinite(g_short)
        
        # Very long time (should give larger g-value)
        t_long = 365 * 24 * 3600  # 1 year
        g_long = fls(t_long, H, r, D, ks, Cs)
        @test g_long > g_short  # g-function increases with time
    end
end
