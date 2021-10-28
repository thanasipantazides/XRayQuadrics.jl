# __precompile__()
module XRayTrace

using Distributed, DistributedArrays
# using BigFloat
# using SharedArrays
using LinearAlgebra
using Random
import JSON, CSV, DataFrames
using Colors
using GLMakie
# using Plots
# using ArbNumerics
import SpecialFunctions

export Particle, PixelatedAttenuator
export Quadric, Plane, Cylinder, Paraboloid, Hyperboloid, TruncatedQuadric, classify, normal, changerepresentation
export inout, batchphotons, transmissionprobability, absorptionprobability, lengthsinattenuator, cylinderentryexit
export PlotTruncatedQuadric, convert_arguments, cartesian_grid, get_mesh, plot, plot!

include("quadrics.jl")
include("optical.jl")
include("tracing.jl")
include("plotting.jl")

end
