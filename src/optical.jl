using XRayQuadrics

struct Particle
    r0::Vector{Float64}     # initial position
    v::Vector{Float64}      # propagation direction (unit)
    E::Float64              # energy
    id::Int64               # id number for tracking
    Particle(r0,v,E,id) = begin
        if length(r0) != 3 || length(v) != 3
            error("3-vectors please")
        else
            new(r0,v/norm(v),E,id)
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
    holes       # a list of cylinders
    largeradius # largest hole radius in attenuator 
    density     # in kg/m^3
    toppoint    # a point (x,y,z) on the one side of the attenuator (in m)
    bottompoint # a point on the other side of the attenuator
    normal      # normal vector to attenuator surface    
    massattenuation # table of mass attenuation coefficients per energy
    bbox        # tuple of two 3-vectors defining corners of rectangular prism bounding the attenuator.
    pitch       # mean hole center-to-hole center spacing
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

function interactiontimes(s::Plane, p::Particle)
    snum = -s.a'*(p.r0 - s.c)
    sdenom = s.a'*p.v
    if sdenom == 0
        return Inf
    else
        return snum/sdenom
    end
end

function interactiontimes(tq::TruncatedQuadric, p::Particle)
    
    rh = [p.r0; 1]
    vh = [p.v; 0]
    Q = tq.q.Q
    qdenom = vh'*Q*vh
    qpnum = -rh'*Q*vh + sqrt(Complex((rh'*Q*vh)^2 - (rh'*Q*rh)*(vh'*Q*vh)))
    qnnum = -rh'*Q*vh - sqrt(Complex((rh'*Q*vh)^2 - (rh'*Q*rh)*(vh'*Q*vh)))
    
    # solution pair is (qpnum/qdenom, qnnum/qdenom)
    
    c1time = interactiontimes(tq.p[1], p)
    c2time = interactiontimes(tq.p[2], p)
    
    if qdenom == 0
        # parallel rays to quadric axis; need to check endcaps for solution.
        if inside(tq.q, c1time*p.v + p.r0) || inside(tq.q, c2time*p.v + p.r0)
            timec = sort!([c1time, c2time])
            return (timec[1], timec[2])
        else
            return (0.0,0.0)
        end
    elseif isinf(c1time) || isinf(c2time) && imag(qpnum) == 0 && imag(qnnum) == 0
        # parallel rays to endcaps; need to check where along quadric it intersects.
        qpnum = real(qpnum)
        qnnum = real(qnnum)
        test_times = sort!([qpnum/qdenom, qnnum/qdenom])
        
        if inside(tq.p[1], tq.p[2], test_times[1]*p.v + p.r0) && inside(tq.p[1], tq.p[2], test_times[2]*p.v + p.r0)
            return (test_times[1], test_times[2])
        else
            return (0.0,0.0)
        end
    end
    
    if imag(qpnum) == 0 && imag(qnnum) == 0
        qpnum = real(qpnum)
        qnnum = real(qnnum)
        timeq = sort!([qpnum/qdenom, qnnum/qdenom])
        timec = sort!([c1time, c2time])
        
        if timec[2] < timeq[1] || timeq[2] < timec[1]
            # hits both planes before quadric, or quadric before both planes
            return (0.0,0.0)
        else
            allt = sort!([timeq;timec])
            return (allt[2], allt[3])
        end
    else
        # no intersection
        return (0.0,0.0)
    end
end

function interactionlength(tq::TruncatedQuadric, p::Particle)
    # how much distance a particle spends inside the TruncatedQuadric
    t1, t2 = interactiontimes(tq, p)
    return norm(t2 - t1)
end