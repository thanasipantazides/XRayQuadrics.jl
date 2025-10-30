using XRayQuadrics
using LinearAlgebra
using GLMakie

function oblique_points(tq::TruncatedQuadric)
    ang = pi/2
    transform = [1 0 0; 0 cos(ang) -sin(ang); 0 sin(ang) cos(ang)]
    dir = [0;0;1]
    v = transform'*dir
    R = 3
    orig = 0.5*(tq.p[1].c .+ tq.p[2].c)  - v*R
    
    n = 50
    points = Point3f[]
    kscale = 0.1
    for k in 1:n
        push!(points, Point3f(orig .+ k*kscale*v))
    end
    return points
end

function intersections()
    # make a cylinder
    a = rand(3)
    # a = [1;1;1]
    # a = [0;0;1]
    a = a/norm(a)
    
    c = rand(3)*2 .- 1
    # c = [5;5;5]
    R = 1
    cyl = Cylinder(
        R,  # radius
        c,
        a
    )
    
    plane1 = Plane(c - 0.5*a, a)
    plane2 = Plane(c + 0.5*a, a)
    
    ccap1 = Point3f(cyl.c + cyl.a'*(plane1.c - cyl.c)*cyl.a)
    ccap2 = Point3f(cyl.c + cyl.a'*(plane2.c - cyl.c)*cyl.a)
    ccap1alt = Point3f(axis_plane_intersection(cyl, plane1))
    ccap2alt = Point3f(axis_plane_intersection(cyl, plane2))
    
    println(ccap1 .- ccap1alt)
    println(ccap2 .- ccap2alt)
    midcyl = Point3f(1/2*(ccap1 .+ ccap2))
    
    tq = TruncatedQuadric(cyl, [plane1,plane2], true)
    
    # make some rays
    n = Int(1e4)
    p = Particle[]
    itimes = zeros(2,n)
    ipos = Point3f[]
    for i in 1:n
        pos = randr()*[0;0;1]*R*3 # random ball of initial positions
        vel = midcyl - pos
        vel = vel/norm(vel)
        thisp = Particle(
            pos,
            vel,
            rand()*10000,
            i
        )
        push!(p, thisp)
        
        itimes[1,i], itimes[2,i] = interactiontimes(tq, thisp)
        if itimes[1,i] != 0 && itimes[2,i] != 0
            push!(
                ipos,
                Point3f(thisp.r0 + thisp.v*itimes[1,i])
            )
            push!(
                ipos,
                Point3f(thisp.r0 + thisp.v*itimes[2,i])
            )
        else
            @warn "missed!"
        end
    end
    
    GLMakie.activate!()
    fig = GLMakie.Figure(size=(500, 800))
    layout3d = GLMakie.GridLayout(fig[1,1])
    scene = GLMakie.LScene(
        layout3d[1,1],
        show_axis=true,
    )
    # plot the surface:
    XRayQuadrics.plot!(scene, tq, caps=false)
    # plot photons:
    # XRayQuadrics.plot!(scene, p, length_scale=0.1)
    # plot some cap and body midpoints
    scatter!(scene, [Point3f(0.0), ccap1, midcyl, ccap2], color=:black, markersize=10)
    # plot the intersection points:
    # scatter!(scene, ipos, color=:blue, markersize=2)
    
    points = oblique_points(tq)
    color = RGBAf[]
    for (k,p) in enumerate(points)
        if inside(tq, p)
            push!(color, RGBAf(0.0,1.0,0.0,1.0))
        else
            push!(color, RGBAf(0.0,0.0,0.0,1.0))
        end
    end
    scatter!(scene, points, color=color, markersize=4)
    
    
    
    display(fig)
    
    # println(cyl.c)
    # println(cyl.a)
    # println(plane1.c)
    # println(plane2.c)
    # println(ccap1, ccap1alt)
    # println(ccap2, ccap2alt)
end