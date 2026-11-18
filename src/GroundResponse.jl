module GroundResponse

using SpecialFunctions
using QuadGK
using LinearAlgebra

# Model abstraction: the AbstractGroundModel type, the concrete model structs and the per-model
# `gmatrix` extension point. Included first so its types are available to the kernels below.
include("model_response.jl")

# Ground response kernels
include("infinite_line_source.jl")          # ils  — Ingersol (1948)
include("infinite_cylindrical_source.jl")   # ics  — Carslaw & Jaeger (1959)
include("finite_line_source.jl")            # fls  — Claesson & Javed (2011)
include("moving_infinite_line_source.jl")   # mils — Pasquier & Lamarche (2022)
include("moving_finite_line_source.jl")     # mfls — Guo et al. (2020)

# Borefield layouts and pairwise geometry
include("borefield.jl")

# Exports
# Abstract type — downstream packages dispatch on this
export AbstractGroundModel
export ILSModel, ICSModel, FLSModel, MILSModel, MFLSModel

# Response arrays of a borefield (boreholes or segments)
export gmatrix

# Kernels
export ils, ics, fls, mils, mfls

# Borefield layouts
export borefield_geometry
export borefield, borefield_rectangle, borefield_line, borefield_circle, borefield_L, borefield_U,
    borefield_open_rectangle

end
