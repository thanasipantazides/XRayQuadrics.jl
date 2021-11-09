# __precompile__()
module XRayTrace

using Distributed, DistributedArrays
# using BigFloat
# using SharedArrays
using LinearAlgebra
using Random
import JSON, HDF5, CSV, DataFrames
using Colors
using GLMakie
# using Plots
# using ArbNumerics
import SpecialFunctions

export Particle, PixelatedAttenuator
export Quadric, Plane, Cylinder, Paraboloid, Hyperboloid, TruncatedQuadric, classify, normal, changerepresentation
export in_out, batch_photons, transmission_probability, absorption_probability, lengthsinattenuator, cylinderentryexit
export PlotTruncatedQuadric, convert_arguments, cartesian_grid, get_mesh, plot, plot!
export bin
export get_reflection_data

include("quadrics.jl")
include("optical.jl")
include("tracing.jl")
include("plotting.jl")
include("analysis.jl")
include("utils.jl")

end
