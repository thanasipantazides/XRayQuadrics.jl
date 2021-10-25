
struct Particle
    r0  # initial position
    v   # propagation direction (unit)
    E   # energy
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
    caps::Vector{Plane}
    surf::Quadric
end

struct WolterIShell
    hyp::Mirror
    par::Mirror
end

struct Optic
    shells::Vector{WolterIShell}
    blocks::Vector{Attenuator}
end

function Particle(r0::Vector{Real}, v::Vector{Real}, E::Real)
    if length(r0) != 3 || length(v) != 3
        error("3-vectors, please")
    end

    return Particle(r0, v/norm(v), E)
end