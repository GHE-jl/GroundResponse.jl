module GroundResponse

using SpecialFunctions
using QuadGK
using PCHIPInterpolation
using LinearAlgebra

# Infinite line source of Ingersol (1948)
include("ground_models/infinite_line_source.jl")
# Infinite cylindrical source of Carslaw and Jaeger (1959)
include("ground_models/infinite_cylindrical_source.jl")
# Finite line source of Claesson and Javed (2011)
include("ground_models/finite_line_source.jl")
# Moving infinite line source of Pasquier and Lamarche (2022)
include("ground_models/moving_infinite_line_source.jl")
# Moving finite line source of Guo et al. (2021)
include("ground_models/moving_finite_line_source.jl")
# ANN-based ground model
include("ground_models/gST_ANN.jl")

# Spatial superposition techniques
include("spatial_superposition.jl")

# Utilities
include("utils.jl")

# Ground model exports
export ils, ics, fls, mils, mfls

# Spatial superposition exports
export bloc_matrix, successive_flux

# Utility exports
export pchip_interpolation, set_nodes, borefield_xy, borefield_radius

end # module GroundResponse
