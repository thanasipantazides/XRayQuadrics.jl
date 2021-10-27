using XRayTrace

"""
    inout(p::Particle, q::Quadric)

Computes entry and exit times for a `Particle` passing through a `Quadric` surface.
"""
function inout(p::Particle, q::Quadric)
    # Project to homogeneous coordinate
    rh = [p.r0;1]   # AN ISSUE WITH THIS: if r = r0 + vΔt, and r0 and v have 1 as their homogenous coord, get (1 + Δt) as the          homogeneous coord for r. FIX
    vh = [p.v;0]    # REMOVE TIME TERM IN HOMEOGENEOUS COORDS?

    Q = q.Q
    # check for coplanar/axial alignment
    denom = vh'*Q*vh
    if denom == 0.0
        if all(Q[1:3, 1:3] .== 0)
            # a plane
            Qd = 2*Q[1:3,end]
            if Qd'*p.v == 0
                return (0,0)
            else
                return (0, (-Qd'*p.r0 - Q[end])/(Qd'*p.v))
            end
        else
            # colinear with nondegenerate quadric axis
            return (0,0)
        end
    end

    post = (sqrt(Complex((vh'*Q*rh)^2 - (vh'*Q*vh)*(rh'*Q*rh))))/denom
    pred = (-vh'*Q*rh)/denom

    if imag(post) == 0
        post = real(post)
        # sort roots in ascending order
        return extrema([pred - post, pred + post])
    else
        # imaginary component implies no intersection
        return (0,0)
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
    absorptionprobability(photon, attenuator)

Return the probability of a photon::Particle being absorbed in attenuator::PixelatedAttenuator.
"""
function absorptionprobability(photon, attenuator)
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
    transmissionprobability(photon, attenuator)

Return the probability of a photon::Particle being transmitted through attenuator::PixelatedAttenuator.
"""
function transmissionprobability(photon, attenuator)
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
    batchphotons(photons, attenuator)

Compute absorption likelihood for a set of `photons` passing through `attenuator`. If multiple processes are available, computation will be parallelized. To execute in parallel, use the `Distributed` module in the calling context, add processes using `addprocs(n)`, and include this source with the `@everywhere` macro before `using` the module:
    
    @everywhere include("Attenuator3D.jl")
    using .Attenuator3D
"""
function batchphotons(photons, attenuator)
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
                localtransmitlikelihood[i] = transmissionprobability(inphotons[i,locali[2][1]], attenuator)
            end
        end

        return Array(transmitlikelihood[1:(end - spares)])

    else
        transmitlikelihood = zeros(size(photons))

        @inbounds for i = 1:length(photons)
            transmitlikelihood[i] = transmissionprobability(photons[i], attenuator)
        end
        return transmitlikelihood
    end
end