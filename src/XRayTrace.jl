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

include("quadrics.jl")
include("optical.jl")
include("tracing.jl")

export Quadric, Plane, Cylinder, Paraboloid, Hyperboloid, inout

export Particle, PixelatedAttenuator
export batchphotons, transmissionprobability, absorptionprobability, lengthsinattenuator, cylinderentryexit

end
