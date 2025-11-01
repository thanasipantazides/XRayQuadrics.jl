using XRayQuadrics
using LinearAlgebra
using GLMakie
import GeometryBasics

function get_intersections(shape::Union{Cylinder, Cone, Paraboloid, Ellipsoid, Hyperboloid}, planes::Vector{Plane}, n::Int)
    
    ca1 = axis_plane_intersection(shape, planes[1])
    ca2 = axis_plane_intersection(shape, planes[2])
    
    mid = (ca1 .+ ca2) / 2
    
    R = 10
    p = Particle[]
    ipos = Point3f[]
    tq = TruncatedQuadric(shape, planes, true)
    for k = 1:n
        pos = randr()*[0;0;1]*R*3
        vel = mid - pos
        vel = vel/norm(vel)
        thisp = Particle(
            pos,
            vel,
            rand()*10000,
            k
        )
        push!(p, thisp)
        itimes1, itimes2 = interactiontimes(tq, thisp)
        if itimes1 != 0 && itimes2 != 0
            push!(
                ipos,
                Point3f(thisp.r0 + thisp.v*itimes1)
            )
            push!(
                ipos,
                Point3f(thisp.r0 + thisp.v*itimes2)
            )
        end
    end
    return ipos
end

function main()
    # using CairoMakie
    GLMakie.activate!(inline=false, title="XRayQuadric shapes")
    # CairoMakie.activate!()

    fig = Figure(size=(1200, 800))

    nax = 3

    center = [1;1;1]
    R = 1
    h = 1.5
    b = 1
    d = 2
    e = 1.5

    for i = 1:nax
        # scene = LScene(fig[i, 1], show_axis=false)
        # ax = [1; cos((i - 1)/nax*pi/2); sin((i - 1)/nax*pi/2)]
        ax = rand(3) .- 0.5
        ax = ax/norm(ax)
        
        center = rand(3)*6 .- 3

        cy = Cylinder(R, center, ax)
        co = Cone(h, center, ax)
        pa = Paraboloid(b, center, ax)
        hy = Hyperboloid(R, b, center, ax)
        el = Ellipsoid(d, e, center, ax)

        p1 = Plane(cy.c - cy.a*0.5, ax)
        p2 = Plane(cy.c + cy.a*1.5, ax)

        surfs = [cy, co, pa, hy, el]

        for j = 1:5
            ipos = get_intersections(surfs[j], [p1, p2], 1000)
            
            q = Quadric(surfs[j])
            tq = TruncatedQuadric(q, [p1, p2], true)
            ax = LScene(fig[i, j], show_axis=false)
            
            scatter!(ax, ipos, color=:blue, markersize=2)
            XRayQuadrics.plot!(ax, tq)
            
            cam3d!(
                ax,
                projectiontype="Orthographic",
                lookat= GeometryBasics.Vec3d(0,0,0),
                upvector=GeometryBasics.Vec3d(0,0,1),
                eyepos= GeometryBasics.Vec3d(10,0,0)
            )
        end
    end

    display(fig)
    # save("test.png", fig; px_per_unit=5)
end