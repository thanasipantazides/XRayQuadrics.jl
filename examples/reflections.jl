using XRayQuadrics
using GLMakie
using LinearAlgebra

function make_shapes()
    center = [0;0;0]
    ax = [1;1;1]
    ax = ax/norm(ax)
    R = 1
    h = 1.5
    b = 1
    d = 2
    e = 1
    
    cy = Cylinder(R, center, ax)
    co = Cone(h, center, ax)
    pa = Paraboloid(b, center, ax)
    hy = Hyperboloid(R, b, center, ax)
    el = Ellipsoid(d, e, center, ax)
    
    quadrics = [cy, co, pa, hy, el]

    p1 = Plane(center .+ ax/2, ax)
    p2 = Plane(center .+ 2*ax, ax)
    
    res = TruncatedQuadric[]
    for q in quadrics
        push!(res, TruncatedQuadric(q, [p1, p2], true))
    end
    return res
end

function make_rays(shape::TruncatedQuadric)
    skew = pi/10
    
end

function parabola_coeffs(r1::Point2d, r2::Point2d)
    # a parabola of the form x = y^2 + ay + b
    a = (r1[1]^2 - r2[1]^2 + r2[2] - r1[2])/(r2[1] - r1[1])
    b = r1[2] - r1[1]^2 - a*r1[1]
    return (a, b)
end

function fit_parabola_coeffs(r::Vector{Point2d})
    A = zeros(length(r),3)
    b = zeros(length(r))
    for k in 1:length(r)
        b[k] = r[k][1] # x-coord of r
        A[k,1] = r[k][2]^2
        A[k,2] = r[k][2]
        A[k,3] = 1
    end
    coeffs = pinv(A)*b
    return coeffs
end

# function fit_hyperbola_coeffs(r::Vector{Point2d})
#     A = zeros(length(r),3)
#     b = zeros(length(r))
#     for k in 1:length(r)
#         b[k] = r[k][1]^2 # x-coord of r
#         A[k,1] = r[k][2]^2
#         A[k,2] = r[k][2]
#         A[k,3] = 1
#     end
#     println(A)
#     println(b)
#     coeffs = pinv(A)*b
    
#     α = coeffs[1]
#     β = coeffs[2]
#     γ = coeffs[3]
    
#     c = β/2/α
#     b2 = α*c^2 - γ
#     a2 = b2/α

#     return [a2,b2,c]
# end

# function fit_hyperbola_coeffs(r::Vector{Point2d})
#     A = zeros(length(r),3)
#     b = zeros(length(r))
#     for k in 1:length(r)
#         b[k] = -r[k][1]^2 # x-coord of r
#         A[k,1] = r[k][1]
#         A[k,2] = r[k][2]^2
#         A[k,3] = 1
#     end
#     println(A)
#     println(b)
#     coeffs = pinv(A)*b
    
#     α = coeffs[1]
#     β = coeffs[2]
#     γ = coeffs[3]
    
#     c = -α/2
#     b2 = (γ - c^2)/β
#     a2 = -β*b2

#     return [a2,b2,c]
# end

function fit_hyperbola_coeffs(r::Vector{Point2d})
    A = zeros(length(r),3)
    b = zeros(length(r))
    for k in 1:length(r)
        b[k] = r[k][2]^2 # x-coord of r
        A[k,1] = r[k][1]^2
        A[k,2] = r[k][1]
        A[k,3] = 1
    end
    println(A)
    println(b)
    coeffs = pinv(A)*b
    
    α = coeffs[1]
    β = coeffs[2]
    γ = coeffs[3]
    
    c = -β/2/α
    b2 = α*c^2 - γ
    a2 = b2/α
    return [a2,b2,c]
end

function hyperbola_coeffs(r1::Point2d, r2::Point2d)
    # a hyperbola of the form x^2/a^2 - y^2/b^2 = 1
    # a2 = a^2, b2 = b^2
    println("denom: ", (r2[1]/r1[1])^2)
    b2 = ((r2[1]/r1[1])^2*r1[2]^2 - r2[2]^2) / (1 - (r2[1]/r1[1])^2)
    a2 = r1[1]^2/(1 + (r1[2]^2/b2))
    return (a2, b2)
end

function reflections()
    tqs = make_shapes()
    
    GLMakie.activate!(title="Reflections")
    fig = Figure(size=(900, 600))
    for (k, tq) in enumerate(tqs)
        ax = LScene(fig[1, k], show_axis=false)
        XRayQuadrics.plot!(ax, tq)
    end
    
    axplane = GLMakie.Axis(fig[2,1:length(tqs)])
    ppoints = [Point2d(300, 53.45), Point2d(0, 51.5), Point2d(300, -53.45), Point2d(0, -51.5)]
    hpoints = [Point2d(0, 51.5), Point2d(-300, 45.6), Point2d(0, -51.5), Point2d(-300, -45.6)]
    
    for k in eachindex(ppoints)
        # ppoints[k] += Point2d(2000,0)
    end
    # for k in eachindex(hpoints)
    #     hpoints[k] += Point2d(2000,0)
    # end
    
    
    # pa,pb = parabola_coeffs(ppoints[1], ppoints[2])
    # ha2,hb2 = hyperbola_coeffs(hpoints[1], hpoints[2])
    
    pcoeffs = fit_parabola_coeffs(ppoints)
    hcoeffs = fit_hyperbola_coeffs(ppoints)
    println(pcoeffs)
    println(hcoeffs)
    
    varyp = -100:0.1:100
    # varxp = varyp.^2 .+ pa.*varyp .+ pb
    varxp = pcoeffs[1].*varyp.^2  .+ pcoeffs[2].*varyp .+ pcoeffs[3]
    varxh = -300:0.1:(0)
    varyh = zeros(length(varxh))
    # for (x,y) in zip(varxh, varyh)
    #     y = sqrt(-hcoeffs[2])*sqrt((x - hcoeffs[3])^2/hcoeffs[1] - 1)
    # end
    
    # VanSpeybroeck + Chase:
    e = 1.0001440190104969
    d = 0.28805876545665093
    f = -2*e^2*d/(e^2 - 1)
    varyh = e^2*(d .+ varxh .+ 2000).^2 .- (varxh .+ 2000).^2
    # varyh = [sqrt(hcoeffs[2])*sqrt((x - hcoeffs[3])^2/hcoeffs[1] - 1) for x in varxh]
    # varyh = sqrt(hcoeffs[2]) .* sqrt.((varxh .- hcoeffs[3]).^2/hcoeffs[1] .- 1)
    # varxh = sqrt(ha2) .* sqrt.(1 .+ (varyp)/hb2)
    
    # println(pa, ", ", pb)
    # println(ha2, ", ", hb2)
    
    lines!(axplane, varxp, varyp, color=:blue, label="parabola")
    lines!(axplane, varxh, varyh, color=:red, label="hyperbola")
    
    pf = Point2d(2000,0)
    scatter!(axplane, [-pf, pf], marker='+', markersize=10, color=:black)
    scatter!(axplane, ppoints, marker='o', markersize=10, color=:blue)
    scatter!(axplane, hpoints, marker='o', markersize=10, color=:red)
    
    display(fig)
end