# __precompile__()
module XRayQuadrics

using Distributed, DistributedArrays
using LinearAlgebra
using Random
import JSON, HDF5, CSV, DataFrames
using Colors
using Makie
# using Plots

export Particle, PixelatedAttenuator
export Quadric, Plane, Cylinder, Paraboloid, Hyperboloid, TruncatedQuadric, normal, changerepresentation
export in_out
export PlotTruncatedQuadric, convert_arguments, cartesian_grid, get_mesh, plot, plot!
export bin
export get_reflection_data
export readsav

include("quadrics.jl")
include("optical.jl")
include("tracing.jl")
# include("plotting.jl")
include("analysis.jl")
include("utils.jl")
include("readsav.jl")

end
