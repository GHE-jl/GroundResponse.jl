"""
    test_moving_infinite_line_source.jl

Tests for the Moving Infinite Line Source (MILS) model.
"""

@testset "Moving Infinite Line Source (MILS)" begin
    # Test parameters
    ks = 3.0                # Ground thermal conductivity [W/mK]
    Cs = 2e6                # Ground volumetric specific heat [J/m³K]
    Cf = 4.2e6              # Groundwater volumetric specific heat [J/m³K]
    r = 0.076               # Radius [m]
    vD = 1e-6               # Darcy velocity [m/s]
    
    @testset "Single time, single radius" begin
        t = 3600.0
        g = mils(t, ks, Cs, Cf, r, vD)
        
        @test g > 0.0
        @test isfinite(g)
    end
    
    @testset "Vector time, single radius" begin
        t = 60:60:3600
        g = mils(t, ks, Cs, Cf, r, vD)
        
        @test isa(g, Vector)
        @test length(g) == length(t)
        @test all(g .> 0.0)
        @test all(isfinite.(g))
        @test all(diff(g) .> 0.0)  # Monotonically increasing
    end
    
    @testset "In-place computation" begin
        t = 60:60:3600
        g_preallocated = Vector{Float64}(undef, length(t))
        
        # Use in-place version
        result = mils!(g_preallocated, t, ks, Cs, Cf, r, vD)
        
        # Result should be the same vector
        @test result === g_preallocated
        
        # Values should be valid
        @test all(g_preallocated .> 0.0)
        @test all(isfinite.(g_preallocated))
        @test all(diff(g_preallocated) .> 0.0)
        
        # Should match the non-inplace version
        g_regular = mils(t, ks, Cs, Cf, r, vD)
        @test all(isapprox.(g_preallocated, g_regular, rtol=1e-10))
    end
    
    @testset "Different Darcy velocities" begin
        t = 3600.0
        
        vD1 = 1e-12  # Very low (impervious)
        vD2 = 1e-6   # Low
        vD3 = 1e-5   # Higher
        
        g1 = mils(t, ks, Cs, Cf, r, vD1)
        g2 = mils(t, ks, Cs, Cf, r, vD2)
        g3 = mils(t, ks, Cs, Cf, r, vD3)
        
        @test all([g1, g2, g3] .> 0.0)
        @test all(isfinite.([g1, g2, g3]))
    end
    
    @testset "Different radii" begin
        t = 3600.0
        
        r1 = 0.05
        r2 = 0.076
        r3 = 0.1
        
        g1 = mils(t, ks, Cs, Cf, r1, vD)
        g2 = mils(t, ks, Cs, Cf, r2, vD)
        g3 = mils(t, ks, Cs, Cf, r3, vD)
        
        @test all([g1, g2, g3] .> 0.0)
        @test all(isfinite.([g1, g2, g3]))
        
        # g-function should decrease with increasing radius
        @test g1 > g2 > g3
    end
    
    @testset "Comparison with ILS" begin
        # With very low velocity, MILS should approach ILS
        t = 3600.0
        vD_low = 1e-12
        
        g_mils = mils(t, ks, Cs, Cf, r, vD_low)
        g_ils = ils(t, r, ks, Cs)
        
        # Both should be valid
        @test g_mils > 0.0
        @test g_ils > 0.0
        
        # They should be different (MILS includes groundwater effects)
        @test g_mils ≠ g_ils
    end
    
    @testset "Type promotion" begin
        t_int = 3600
        g1 = mils(t_int, ks, Cs, Cf, r, vD)
        @test isa(g1, AbstractFloat)
    end
    
    @testset "In-place type consistency" begin
        t = [60.0, 120.0, 180.0]
        g = Vector{Float64}(undef, length(t))
        
        mils!(g, t, ks, Cs, Cf, r, vD)
        
        @test all(isa.(g, Float64))
    end
end
