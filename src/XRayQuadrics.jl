# __precompile__()
module XRayQuadrics

using Distributed, DistributedArrays
using LinearAlgebra, SparseArrays
using Random
import JSON, HDF5, CSV, DataFrames
using Colors
using GLMakie

include("utils.jl")
include("quadrics.jl")
include("optical.jl")
include("tracing.jl")
include("plotting.jl")
include("analysis.jl")

export r_euler3, r_euler2, r_euler1, r_min_arc
export Particle, PixelatedAttenuator
export Quadric, Plane, Cylinder, Cone, Ellipsoid, Paraboloid, Hyperboloid, TruncatedQuadric, normal, changerepresentation
export in_out, solve_quadratic, batch_photons, batch_photons_through_attenuator
export convert_arguments, cartesian_grid, get_mesh, plot, plot!, axis_plane_intersection, interactionlength, interactiontimes, inside
export Parabola, Hyperbola
export fitp, fith, focus, randr
export bin
export get_reflection_data, get_mass_attenuation_data, get_photon_data, interpolateattenuation, get_empirical_attenuation_data


end
