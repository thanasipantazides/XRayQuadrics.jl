using XRayQuadrics

struct Particle
    r0::Vector{Float64}     # initial position
    v::Vector{Float64}      # propagation direction (unit)
    E::Float64              # energy
    id::Int64               # id number for tracking
    Particle(r0,v,E) = begin
        if length(r0) != 3 || length(v) != 3
            error("3-vectors please")
        else
            new(r0,v/norm(v),E)
        end
    end
end

struct Attenuator   # UNUSED
    density
    toppoint
    bottompoint
    normal
    massattenuation
end

struct PixelatedAttenuator
    holes
    largeradius # largest hole radius in attenuator 
    density     # in kg/m^3
    toppoint    # a point (x,y,z) on the one side of the attenuator (in m)
    bottompoint # a point on the other side of the attenuator
    normal      # normal vector to attenuator surface    
    massattenuation # table of mass attenuation coefficients per energy
end

struct Mirror
    tq::TruncatedQuadric
    reflection_probability::Matrix{Float64}
    matl::String
end

struct WolterIShell
    hyp::Mirror
    par::Mirror
end

struct Optic
    shells::Vector{WolterIShell}
    blocks::Vector{Attenuator}
end

# function Particle(r0::Vector{Real}, v::Vector{Real}, E::Real)
#     if length(r0) != 3 || length(v) != 3
#         error("3-vectors, please")
#     end

#     return Particle(r0, v/norm(v), E)
# end