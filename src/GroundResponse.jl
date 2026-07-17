module GroundResponse

using SpecialFunctions
using QuadGK
using LinearAlgebra
using DSP
using PCHIPInterpolation: Interpolator

# Model abstraction: the AbstractGroundModel type, the concrete model structs, the per-model
# `_borehole_response` extension point, and the high-level `ground_response` dispatcher. Included
# first so its types are available to the kernels and the superposition backends below.
include("model_response.jl")

# Ground response kernels
include("infinite_line_source.jl")          # ils  — Ingersol (1948)
include("infinite_cylindrical_source.jl")   # ics  — Carslaw & Jaeger (1959)
include("finite_line_source.jl")            # fls  — Claesson & Javed (2011)
include("moving_infinite_line_source.jl")   # mils — Pasquier & Lamarche (2022)
include("moving_finite_line_source.jl")     # mfls — Guo et al. (2020)
include("short_term_ann.jl")                # short-term ANN — Pasquier et al. (2018)

# Time sub-sampling / interpolation, then the spatial-superposition backends and borefield layouts.
include("time_sampling.jl")
include("spatial_superposition.jl")
include("borefield.jl")

# Exports
# Abstract type — downstream packages dispatch on this
export AbstractGroundModel
export ILSModel, ICSModel, FLSModel, MILSModel, MFLSModel

# High-level interface
export borehole_response, ground_response

# Backends — available for direct use
export ils, ics, fls, mils, mfls, short_term_response

# Spatial superposition
export bloc_matrix, successive_flux, uniform_flux, segment_response, segment_response_marching

# Borefield layouts
export borefield_geometry
export borefield, borefield_rectangle, borefield_line, borefield_circle, borefield_L, borefield_U,
    borefield_open_rectangle

end
