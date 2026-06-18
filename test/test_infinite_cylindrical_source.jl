"""
    test_infinite_cylindrical_source.jl

Tests for the Infinite Cylindrical Source (ICS) model.
"""

@testset "Infinite Cylindrical Source (ICS)" begin
    # Test parameters
    r = 0.076               # Measurement radius [m]
    rc = 0.076              # Cylinder radius [m]
    ks = 3.0                # Ground thermal conductivity [W/mK]
    Cs = 2e6                # Ground volumetric specific heat [J/m³K]
    
    @testset "Single time, single radius" begin
        t = 3600.0
        g = ics(t, r, rc, ks, Cs)
        
        @test g > 0.0
        @test isfinite(g)
    end
    
    @testset "Vector time, single radius" begin
        t = 60:60:3600
        g = ics(t, r, rc, ks, Cs)
        
        @test isa(g, Vector)
        @test length(g) == length(t)
        @test all(g .> 0.0)
        @test all(isfinite.(g))
        @test all(diff(g) .> 0.0)  # Monotonically increasing
    end
    
    @testset "Single time, vector radius" begin
        t = 3600.0
        r_vec = [0.076, 0.1, 0.15, 0.2]
        g = ics(t, r_vec, rc, ks, Cs)
        
        @test isa(g, Vector)
        @test length(g) == length(r_vec)
        @test all(g .> 0.0)
        @test all(isfinite.(g))
        
        # g-function should decrease with increasing radius
        @test all(diff(g) .< 0.0)
    end
    
    @testset "Vector time, vector radius (2D matrix)" begin
        t = 60:60:3600
        r_vec = [0.076, 0.1, 0.15]
        g = ics(t, r_vec, rc, ks, Cs)
        
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
    
    @testset "Different cylinder radii" begin
        t = 3600.0
        r_meas = 0.1
        
        # Test with different cylinder radii
        rc1 = 0.05
        rc2 = 0.10
        rc3 = 0.15
        
        g1 = ics(t, r_meas, rc1, ks, Cs)
        g2 = ics(t, r_meas, rc2, ks, Cs)
        g3 = ics(t, r_meas, rc3, ks, Cs)
        
        # All should be valid g-functions
        @test all([g1, g2, g3] .> 0.0)
        @test all(isfinite.([g1, g2, g3]))
    end
    
    @testset "Type promotion" begin
        t_int = 3600
        g1 = ics(t_int, r, rc, ks, Cs)
        @test isa(g1, AbstractFloat)
    end
end
