using XRayQuadrics

"""
    in_out(p::Particle, q::Quadric)

Computes entry and exit times for a `Particle` passing through a `Quadric` surface. Returns an Array of length 2, with the entry and exit times in increasing order, unless there is only one intersection with the surface, in which case the time is the first element. If there is no interaction, returns two zeros.
"""
function in_out(p::Particle, q::Quadric)
    # Project to homogeneous coordinate
    rh = [p.r0;1]   # AN ISSUE WITH THIS: if r = r0 + vΔt, and r0 and v have 1 as their     homogenous coord, get (1 + Δt) as the homogeneous coord for r. FIX
    vh = [p.v;0]    # REMOVE TIME TERM IN HOMEOGENEOUS COORDS?

    Q = q.Q
    Qr = Q[1:3,1:3]
    qd = Q[1:3,end]
    q0 = Q[end]

    # check for coplanar/axial alignment
    denom = vh'*Q*vh

    if all(Qr .== 0)                # plane
        if qd'*p.v == 0             # ray parallel to plane
            return [0, 0]
        else                        # ray intersects plane
            a = 2*qd
            return [(-q0 - a'*p.r0)/(a'*p.v), 0]
        end
    else                            # nonreducible quadric
        if denom == 0               # ray parallel to axis, not hyperboloid
            if qd'*p.v != 0         # paraboloid
                return [(p.r0'*Qr*p.r0 + 2*qd'*p.r0 + q0)/(2*qd'*p.v), 0]
            else                    # cylinder, but ray parallel to axis
                return [0, 0]
            end
        else                        # ray oblique to axis, or hyperboloid
            post = (sqrt(Complex((vh'*Q*rh)^2 - (vh'*Q*vh)*(rh'*Q*rh))))/denom
            pred = (-vh'*Q*rh)/denom
            if imag(post) == 0      # ray hits quadric
                post = real(post)
                return [min(pred - post, pred + post), max(pred - post, pred + post)]
            else                    # imaginary solution ⟹ no intersection with quadric
                return [0, 0]
            end
        end
    end

    # if denom == 0.0
    #     if all(Q[1:3, 1:3] .== 0)
    #         # a plane
    #         Qd = 2*Q[1:3,end]
    #         if Qd'*p.v == 0
    #             return (0,0)
    #         else
    #             return (0, (-Qd'*p.r0 - Q[end])/(Qd'*p.v))
    #         end
    #     else
    #         # colinear with nondegenerate quadric axis
    #         return (0,0)
    #     end
    # end

    # post = (sqrt(Complex((vh'*Q*rh)^2 - (vh'*Q*vh)*(rh'*Q*rh))))/denom
    # pred = (-vh'*Q*rh)/denom

    # if imag(post) == 0
    #     post = real(post)
    #     # sort roots in ascending order
    #     return extrema([pred - post, pred + post])
    # else
    #     # imaginary component implies no intersection
    #     return (0,0)
    # end
end

function in_out_diagonalize(c::Float64, pblock::Matrix, Vblock::Matrix, Qblock::Matrix)
    
end

"""
    intersect_photons_with_quadrics(p::Array{Particle}, tq::TruncatedQuadric) 

Computes intersection positions of a set of photons with a set of surfaces.
"""
function intersect_photons_with_quadrics(p::Vector{Particle}, tq::Vector{TruncatedQuadric})
    n_surfaces = length(tq)
    n_photons = length(p)

    # matrix structure: grouped by surface (not grouped by photon)

    γ = zeros(n_photons*n_surfaces)
    q = zeros(n_photons*n_surfaces, n_photons*n_surfaces)
    for i in 1:length(q[1,:])
        # γ[i] = p[i]
        for j in 1:length(q[:,1])

        end
    end
end

"""
    build_tracing_geometry(p::Vector{Particle}, tq::Vector{TruncatedQuadric})

Builds generic matrix-vector quadratic terms `A`, `b`, `c` from provided `Quadric` surfaces and ray data.
"""
function build_tracing_geometry(p::Vector{Particle}, tq::Vector{TruncatedQuadric})

end

"""
    solve_quadratic(A::Array, b::Vector, c::Real)

Solves quadratic matrix-vector equation. Returns two solution vectors `x` to `x'*A*x + b'*x + c = 0`.
"""
function solve_quadratic(A::Array, b::Vector, c::Real)
    # check `A` symmetric
    is_symmetric = A == A'

    # check size of `b`
    if length(b) != length(A[:,1])
        error("matrix and vector must have same length")
    end

    # check if `A` is invertible:
    use_pseudo_inverse = rank(A) != length(A[:,1])

    # complete the square parameters:
    Ainv = zeros(length(A[:,1]), length(A[1,:]))
    if use_pseudo_inverse
        Ainv = pinv(A)
    else
        Ainv = inv(A)
    end
    h = -0.5*Ainv*b
    k = c - 0.25*b'*Ainv*b

    # diagonalize A:
    (Λ, V) = eigen(A)

    # get solutions in principal axes space
    zp = sqrt.(-k ./ Λ)
    zn = -zp

    # transform solutions back to original space
    xp = inv(V)*zp + h
    xn = inv(V)*zn + h

    return (xp, xn)
end

"""
    interacts(p::Particle, tq::TruncatedQuadric)

!!! warning
    FINISH IMPLEMENTING: should return a reflected/transmitted copy of `p` and a probability (of reflection or transmission, depending on whether `tq` is a mirror or not) 
"""
function interacts(p::Particle, tq::TruncatedQuadric)
    quadric = tq.q
    Qr = quadric.Q[1:3,1:3]
    qd = quadric.Q[1:3,end]
    q0 = quadric.Q[end]

    plane1 = tq.p[1]
    plane2 = tq.p[2]

    qtimes = in_out(p,quadric)
    p1time = in_out(p,plane1)[1]
    p2time = in_out(p,plane2)[1]

    nzeros = length(findall(qtimes .== 0))

    if nzeros == 2                  # parallel to quadric axis
        # project into both planes, check if within quadric at those positions
        # PROBLEM: HOW TO TELL IF IT'S *WITHIN* QUADRIC?

    elseif nzeros == 1              # quadric must be paraboloid
        # badorders:
        #   spp (anti-axis)
        #   pps (with axis)
        #   ppss    CANNOT APPEAR HERE
        #   sspp    CANNOT APPEAR HERE
        if tq.caps
            dir = sign(qd'*p.v)
            if dir == 1             # path parallel to axis

            elseif dir == 0
                error("Particle was supposed to travel parallel to Paraboloid axis.")
            else                    # path antiparallel to axis

            end
        else
            # psp is only passing case
        end
    else                            # hits walls of quadric twice
        
    end
end

"""
    cylinderentryexit(particle, cylinder)

Compute times at which a ray passes through a cylinder. The ray has initial position `particle.r0` and velocity `particle.v`. The cylinder has radius `cylinder.R`, height `cylinder.h`, a point one one of its caps `cylinder.c`, and a unit-length axis pointing from `cylinder.c` towards its other cap `cylinder.a`. 

Two scalars are returned: the time at which the ray enters the cylinder and the time at which it leaves the cylinder. If the ray misses the cylinder, zeros are returned for both entry and exit time.
"""
function cylinderentryexit(particle, cylinder)
    r0 = particle.r0
    v = particle.v
    c = cylinder.c
    a = cylinder.a
    R = cylinder.R
    h = cylinder.h

    ttop = (c'*a - r0'*a)/(v'*a)
    tbottom = (c'*a + h - r0'*a)/(v'*a)

    topcapradius = norm(r0 + v*ttop - c)
    bottomcapradius = norm(r0 + v*tbottom - c - a*h)

    qa = v'*v - (v'*a)^2
    qb = 2*((r0'*v - c'*v) - ((v'*a)*(r0'*a - c'*a)))
    qc = norm(r0 - c)^2 - (r0'*a - c'*a)^2 - R^2

    # particle coaxial with cylinder:
    if qa == 0
        # check if it hits caps, if so return caps times
        if topcapradius ≤ R
            if bottomcapradius ≤ R
                # rectify to zero if particle starts within cylinder
                ttop = max(ttop, 0)
                tbottom = max(tbottom, 0)
                return (min(ttop, tbottom), max(ttop, tbottom))
            else
                error("particle supposedly parallel to cylinder axis")
            end
        else
            # outside cap radius, miss
            if bottomcapradius ≤ R
                error("particle supposedly parallel to cylinder axis")
            end

            return (0, 0)
        end
        
    # photon oblique to cylinder axis:
    else
        root1 = (-qb + sqrt(Complex(qb^2 - 4*qa*qc)))/(2*qa)
        root2 = (-qb - sqrt(Complex(qb^2 - 4*qa*qc)))/(2*qa)

        if imag(root1) == 0
            root1 = real(root1)
            root2 = real(root2)
            # hits infinite cylinder somewhere
            t1 = min(root1, root2)
            t2 = max(root1, root2)

            # list of all intersection times with cap planes and cylinder
            dt = [tbottom, ttop, t1, t2]
            # order intersection times
            I = sortperm(dt)

            # check if first two intersection times are with cylinder (t1, t2)
            if length(intersect(I[1:2], [3, 4])) == 2
                # miss
                return (0, 0)
            end
            # check if last two intersection times are with cylinder (t1, t2)
            if length(intersect(I[3:4], [3, 4])) == 2
                # miss
                return (0, 0)
            end

            # if you've made it this far, cylinder is hit somewhere
            tin = dt[I[2]]
            tout = dt[I[3]]

            tin = max(tin, 0)
            tout = max(tout, 0)
            return (tin, tout)

        # miss altogether (complex roots for intersection with cylinder) 
        else
            return (0, 0)
        end
    end
end

"""
    lengthsinattenuator(photon, attenuator)

Get distances a `photon::Particle` travels within an `attenuator::PixelatedAttenuator` which has cylinders cut out of it. 
"""
function lengthsinattenuator(photon, attenuator)
    topcylinderstimes = Array{Float64}(undef, 0, 2)
    # filter cylinders along path
    hits = falses(length(attenuator.holes[:]))
    for c = 1:length(attenuator.holes[:])
        if norm(cross(photon.v,attenuator.holes[c].a)) == 0
            # if cylinder parallel to ray
            hits[c] = norm(attenuator.holes[c].c - photon.r0 - photon.v*((photon.v'*(attenuator.holes[c].c - photon.r0))./(photon.v'*photon.v))) ≤ attenuator.holes[c].R
        else
            # if cylinder oblique to ray
            hits[c] = cross(photon.v, attenuator.holes[c].a)'*(photon.r0 - attenuator.holes[c].c)/norm(cross(photon.v, attenuator.holes[c].a)) ≤ attenuator.holes[c].R
        end
        
        # if this photon flies within this cylinder's radius, check for collision
        if hits[c]
            (tin, tout) = cylinderentryexit(photon, attenuator.holes[c])

            # make sure it actually hits cylinder
            if tout - tin > 0
                topcylinderstimes = cat(topcylinderstimes, [tin tout], dims=1)
            end
        end
    end
    
    # sort cylinders visited in order of visit time (increasing distance from r0 along v)
    sortedcylinderstimes = zeros(size(topcylinderstimes))
    if length(topcylinderstimes[:,1]) > 1
        topcylinderorder = sortperm(topcylinderstimes[:,1])
        sortedcylinderstimes = topcylinderstimes[topcylinderorder,:]
    else
        sortedcylinderstimes = topcylinderstimes
    end
    
    # times at which photon hits attenuator top surface and bottom surface
    tslabbottom = (attenuator.bottompoint'*attenuator.normal - photon.r0'*attenuator.normal)/(photon.v'*attenuator.normal)
    tslabtop = (attenuator.toppoint'*attenuator.normal - photon.r0'*attenuator.normal)/(photon.v'*attenuator.normal)

    tslabin = min(tslabbottom, tslabtop)
    tslabout = max(tslabbottom, tslabtop)

    if length(sortedcylinderstimes) > 0
        # times spent within attenuator material (complement of time spent in cylinders)
        timediffs = [sortedcylinderstimes[:,1]; tslabout] - [tslabin; sortedcylinderstimes[:,2]]
    else
        timediffs = tslabout - tslabin
    end

    # path lengths in material
    lengths = norm(photon.v).*timediffs

    return lengths
end

"""
    absorption_probability(photon, attenuator)

Return the probability of a photon::Particle being absorbed in attenuator::PixelatedAttenuator.
"""
function absorption_probability(photon, attenuator)
    lengths = lengthsinattenuator(photon, attenuator)
    if length(lengths) != 0
        transmissionlikelihoods = zeros(BigFloat, length(lengths))
        for i = 1:length(lengths)
            # interpolate attenuation coefficients from table data:
            attenuationcoeff = interpolateattenuation(photon.E, attenuator.massattenuation)
            
            # probability of passing through material:
            transmissionlikelihoods[i] = exp(BigFloat(-lengths[i]*attenuationcoeff))
        end

        # complement of total transmission likelihood is absorption likelihood
        absorptionlikelihood = 1 - prod(transmissionlikelihoods)
        
        return absorptionlikelihood
    else
        return 0.0
    end
end

"""
    transmission_probability(photon, attenuator)

Return the probability of a photon::Particle being transmitted through attenuator::PixelatedAttenuator.
"""
function transmission_probability(photon, attenuator)
    lengths = lengthsinattenuator(photon, attenuator)
    if length(lengths) != 0
        transmissionlikelihoods = zeros(BigFloat, length(lengths))
        for i = 1:length(lengths)
            # interpolate attenuation coefficients from table data:
            attenuationcoeff = interpolateattenuation(photon.E, attenuator.massattenuation)
            
            # probability of passing through material:
            transmissionlikelihoods[i] = exp(BigFloat(-lengths[i]*attenuationcoeff))
        end

        # get total transmission likelihood
        return prod(transmissionlikelihoods)
        
    else
        return 1.0
    end
end

"""
    batch_photons(photons, attenuator)

Compute absorption likelihood for a set of `photons` passing through `attenuator`. If multiple processes are available, computation will be parallelized. To execute in parallel, use the `Distributed` module in the calling context, add processes using `addprocs(n)`, and include this source with the `@everywhere` macro before `using` the module:
    
    @everywhere include("Attenuator3D.jl")
    using .Attenuator3D
"""
function batch_photons(photons, attenuator)
    # check if multiple processes available:
    if nprocs() > 1
        # need to factor jobsize by nworkers()
        spares = nworkers() - rem(length(photons), nworkers())
        # populate extra "dummy" photons to make array correct size
        appendix = Array{Particle}(undef, spares, 1)
        fill!(appendix, photons[1])
        inphotons = [photons[:]; appendix]
        
        # size of batch to run on each worker
        batchsize = length(inphotons[:]) ÷ nworkers()

        if rem(length(inphotons[:]), nworkers()) != 0
            error("Cannot distribute job to number of workers")
        end

        # shape array so each worker gets one column
        inphotons = reshape(inphotons, batchsize, nworkers())
        # make distributed array for result on all workers
        transmitlikelihood = dzeros(BigFloat, size(inphotons), workers(), [1, nworkers()])

        @sync @distributed for j in 1:nworkers()
            localtransmitlikelihood = localpart(transmitlikelihood)
            locali = DistributedArrays.localindices(transmitlikelihood)

            for i = 1:length(localtransmitlikelihood)
                localtransmitlikelihood[i] = transmission_probability(inphotons[i,locali[2][1]], attenuator)
            end
        end

        return Array(transmitlikelihood[1:(end - spares)])

    else
        transmitlikelihood = zeros(size(photons))

        @inbounds for i = 1:length(photons)
            transmitlikelihood[i] = transmission_probability(photons[i], attenuator)
        end
        return transmitlikelihood
    end
end