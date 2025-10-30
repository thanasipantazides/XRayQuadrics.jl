using XRayQuadrics
using LinearAlgebra
using GLMakie
import GeometryBasics

function main()
    # using CairoMakie
    GLMakie.activate!(inline=false)
    # CairoMakie.activate!()

    fig = Figure(size=(1200, 800))

    nax = 5

    center = [1;1;1]
    R = 1
    h = 1.5
    b = 1
    d = 2
    e = 1

    for i = 1:nax
        # scene = LScene(fig[i, 1], show_axis=false)
        # ax = [1; cos((i - 1)/nax*pi/2); sin((i - 1)/nax*pi/2)]
        ax = [1; 0; sin((i-1)/nax*pi/2)]
        ax = ax/norm(ax)

        cy = Cylinder(R, center, ax)
        co = Cone(h, center, ax)
        pa = Paraboloid(b, center, ax)
        hy = Hyperboloid(R, b, center, ax)
        el = Ellipsoid(d, e, center, ax)

        p1 = Plane(cy.a, ax)
        p2 = Plane(cy.a*4, ax)

        surfs = [cy, co, pa, hy, el]

        for j = 1:5
            q = Quadric(surfs[j])
            tq = TruncatedQuadric(q, [p1, p2], true)
            ax = LScene(fig[i, j], show_axis=false)
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

    display(GLMakie.Screen(), fig)
    # save("test.png", fig; px_per_unit=5)
end