__precompile__()
module XRayTrace

using Distributed, DistributedArrays
# using BigFloat
# using SharedArrays
using LinearAlgebra
using Random
import JSON, CSV, DataFrames
using Plots
# using ArbNumerics
import SpecialFunctions

export Particle, PixelatedAttenuator
export Quadric, Plane, Cylinder, Paraboloid, Hyperboloid, classify, changerepresentation
export inout, batchphotons, transmissionprobability, absorptionprobability, lengthsinattenuator, cylinderentryexit

include("quadrics.jl")
include("optical.jl")
include("tracing.jl")

end
