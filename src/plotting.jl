using XRayQuadrics
# import Makie: plot, plot!
using GLMakie
using Colors

# macro expands to:
#   const PlotTQ{ArgTypes} = Combined{plottq, ArgTypes}
#   plottq(args...; kw_args...) = ...
#   plottq!(args...; kw_args...) = ...
#   function default_theme(scene, ::PlotTQ) = ...

# @Makie.recipe(PlotTQ, obj) do scene
#     Makie.Attributes(;
#         plot_caps = true,
#         color = :red,
#         seriesalpha = 0.667,
#         light_source = [0;0;1]
#     )
# end

# Makie.plottype(::TruncatedQuadric) = PlotTQ{<:Tuple{TruncatedQuadric}}

# function Makie.plot!(plot::PlotTQ{<:Tuple{TruncatedQuadric}})
#     tq = plot[:obj][]

#     (Xs, Ys, Zs) = cartesian_grid(tq)

#     X = Xs[1]
#     Y = Ys[1]
#     Z = Zs[1]

#     Makie.surface!(
#         plot, X, Y, Z,
#         shading=false,
#         color=fill(RGBA(1.0,0.0,0.0,0.4), size(X)...)
#     )
#     if plot[:plot_caps][]
#         Makie.surface!(
#             plot, Xs[2], Ys[2], Zs[2],
#             shading=false,
#             color=fill(RGBA(0.0,0.0,1.0,0.4), size(Xs[2]))
#         )
#         Makie.surface!(
#             plot, Xs[3], Ys[3], Zs[3],
#             shading=false,
#             color=fill(RGBA(0.0,0.0,1.0,0.4), size(Xs[3]))
#         )
#     end
#     return plot
# end

# Plots.@recipe function f(tq::TruncatedQuadric; plot_caps = false)
#     # color       --> :blue
#     seriesalpha --> 0.667
#     seriestype  :=  :surface
#     legend      --> false
#     # markershape --> (add_marker ? :circle : :none)
#     delete!(plotattributes, :add_marker)
#     (X, Y, Z) = cartesian_grid(tq)
#     @series begin
#         seriestype  := :surface
#         legend      --> false
#         X[1], Y[1], Z[1]
#     end
#     # if plot_caps
#         @series begin
#             seriestype  := :surface
#             legend      --> false
#             X[2], Y[2], Z[2]
#         end
#         @series begin
#             seriestype  := :surface
#             legend      --> false
#             X[3], Y[3], Z[3]
#         end
#     # end
# end

# @Makie.recipe(TruncatedQuadricPlot, tq) do scene
#     Theme(;
#         color = :red,
#     )
# end

# function Makie.plot!(tq_plot::TruncatedQuadricPlot{<:Tuple{TruncatedQuadric}})
#     println("Makie.plot!(tq_data::TruncatedQuadricPlot) ran.")

#     (X, Y, Z) = cartesian_grid(tq_plot[:tq])

#     surface!(tq_plot, X, Y, Z)

#     return tq_plot
# end

# function plot(ax::Makie.LScene, tq::TruncatedQuadric)
#     (X, Y, Z) = cartesian_grid(tq)

#     surface!(ax, X, Y, Z)
# end

function plot!(ax::Makie.LScene, ps::Vector{Particle}; length_scale=1e-6)
    n = length(ps)
    
    us = Point3f[(0.0,0.0,0.0) for k in 1:3*n]
    
    k = 1
    korig = 1
    while k < length(us)
        us[k] = ps[korig].r0
        us[k+1] = ps[korig].r0 .+ length_scale*ps[korig].v
        us[k+2] = Point3f(NaN)
        
        korig += 1
        k += 3
    end
    
    lines!(ax, us, color=:blue, linewidth=0.5)
end

function plot!(ax::GLMakie.LScene, tq::TruncatedQuadric; caps=false)
    (X, Y, Z, T, debug) = cartesian_grid(tq)

    colors = Dict(Paraboloid=> :green, Hyperboloid=> :blue, Cylinder=> :red, Cone=> :yellow, Ellipsoid=> :magenta)

    caps = tq.caps
    surface!(
        ax, 
        X[1], Y[1], Z[1],
        colorrange=(-30, -20),
        highclip=(colors[T], 0.5),
        transparency=true
    )
    
    if caps
        # lines!(ax, debug["origin"], linewidth=5, color=:black)
        surface!(
            ax, 
            X[2], Y[2], Z[2],
            # color=colors[T]
            colorrange=(-30, -20),
            highclip=(colors[T], 0.3),
            transparency=true
        )
        surface!(
            ax, 
            X[3], Y[3], Z[3],
            # color=colors[T]
            colorrange=(-30, -20),
            highclip=(colors[T], 0.3),
            transparency=true
        )
    end
end


function cartesian_grid(tq::TruncatedQuadric)
    q = tq.q
    ps = tq.p
    s = changerepresentation(q)
    T = typeof(s)

    (nθ, nζ) = (30, 100)
    # (Xs, Ys, Zs) = get_mesh(s, ps, nθ, nζ)
    (Xs, Ys, Zs, debug) = get_mesh(tq, nθ, nζ)
    return (Xs, Ys, Zs, T, debug)
end

function interactionsites(tq::TruncatedQuadric, p::Particle)
    t1, t2 = interactiontimes(tq, p)
    if t1 == t2 == 0.0
        return (Point3f(NaN), Point3f(NaN))
    else
        return (Point3f(p.r0 + p.v*t1), Point3f(p.r0 + p.v*t2))
    end
end

function get_mesh(tq::TruncatedQuadric, nθ::Int64, nζ::Int64)
    ca1, ca2 = axis_plane_intersection(tq)
    ax = ca2 - ca1
    ax = ax/norm(ax)
    mid = 0.5*(ca1 + ca2)
    dist = norm(ca2 - ca1)
    ca1proj = ax'*ca1
    ca2proj = ax'*ca2
    
    xordinate = ax - rand(3) # a reference direction normal to the axis
    xordinate = cross(ax, xordinate)
    xordinate = xordinate/norm(xordinate)
    yordinate = cross(ax, xordinate)
    yordinate = yordinate/norm(yordinate)
    outer_transform = [xordinate yordinate ax]
    # outer_transform = r_min_arc(ax, ordinate)
    
    θ = range(0, stop=2*π, length=nθ) # each ray should hit surface twice, so just go up to π
    # ζ = range(ca1proj, stop=ca2proj, length=nζ)
    ζ = range(0.01, stop=0.99, length=nζ)
    
    X = zeros(nθ, nζ)
    Y = zeros(nθ, nζ)
    Z = zeros(nθ, nζ)
    
    debug = Dict("origin"=>Point3f[])
    for j in 1:nζ
        origin = ca1 + ζ[j]*(ca2 - ca1) #+ xordinate*dist*1e-1 # walk along the axis
        push!(debug["origin"], Point3f(origin))
        for i in 1:nθ
            dir = outer_transform*r_euler3(θ[i])*[1;0;0]
            ray = Particle(origin, dir, 0, 0)
            p1, p2 = interactionsites(tq, ray)
            X[i,j] = p1[1]
            Y[i,j] = p1[2]
            Z[i,j] = p1[3]
            # X[2*i,j] = p2[1]
            # Y[2*i,j] = p2[2]
            # Z[2*i,j] = p2[3]
        end
    end
    
    Xc1 = [ca1[1]*ones(nθ, 1) X[:,1]]
    Xc2 = [ca2[1]*ones(nθ, 1) X[:,end]]
    Yc1 = [ca1[2]*ones(nθ, 1) Y[:,1]]
    Yc2 = [ca2[2]*ones(nθ, 1) Y[:,end]]
    Zc1 = [ca1[3]*ones(nθ, 1) Z[:,1]]
    Zc2 = [ca2[3]*ones(nθ, 1) Z[:,end]]
    
    return ((X,Xc1,Xc2), (Y,Yc1,Yc2), (Z,Zc1,Zc2), debug)
end

function get_mesh(s::Union{Cylinder, Paraboloid, Hyperboloid, Cone, Ellipsoid}, ps::Vector{Plane}, nθ::Int64, nζ::Int64)

    ca1 = axis_plane_intersection(s, ps[1])
    ca2 = axis_plane_intersection(s, ps[2])
    
    projaca1 = s.a'*(ca1 - s.c)
    projaca2 = s.a'*(ca2 - s.c)
    # if projaca2 < projaca1
    #     temp = ca1
    #     ca1 = ca2
    #     ca2 = temp
        
    #     # temp = projaca1
    #     # projaca1 = projaca2
    #     # projaca2 = temp
    # end

    h = norm(ca2 - ca1)

    # h1 = ca1'*s.a
    # h2 = ca2'*s.a
    h1 = projaca1
    h2 = projaca2
    # h1 = ca1'*s.a
    # h2 = ca2'*s.a
     
    h1proj = abs(h1)
    h2proj = h1proj + h

    a = s.a
    c = s.c
    θ = range(0, stop=2π, length=nθ)
    # ζ = range(0, stop=h, length=nζ)
    ζ = range(h1proj, stop=h2proj, length=nζ)

    # todo: bound ζ so that it conforms to shape boundary (nonimaginary)

    X = zeros(nθ, nζ)
    Y = zeros(nθ, nζ)
    Z = zeros(nθ, nζ)

    if typeof(s) == Cylinder
        R = s.R
        for i in 1:nθ
            for j in 1:nζ
                X[i,j] = R*cos(θ[i])
                Y[i,j] = R*sin(θ[i])
                Z[i,j] = ζ[j]
            end
        end
    elseif typeof(s) == Cone
        for i in 1:nθ
            for j in 1:nζ
                X[i,j] = 1/s.h*ζ[j]*cos(θ[i])
                Y[i,j] = 1/s.h*ζ[j]*sin(θ[i])
                Z[i,j] = ζ[j]
            end
        end
    elseif typeof(s) == Ellipsoid
        d = s.d
        e = s.e
        for i in 1:nθ
            for j in 1:nζ
                # X[i,j] = d*cos(θ[i])*sin(acos(ζ[j]/h2))
                # Y[i,j] = d*sin(θ[i])*sin(acos(ζ[j]/h2))
                X[i,j] = d*cos(θ[i])*sin(acos((ζ[j] - h1proj)/(h2proj - h1proj)))
                Y[i,j] = d*sin(θ[i])*sin(acos((ζ[j] - h1proj)/(h2proj - h1proj)))
                
                Z[i,j] = ζ[j]
            end
        end
    elseif typeof(s) == Paraboloid
        println("derived center: ", s.c)
        b = s.b
        for i in 1:nθ
            for j in 1:nζ
                X[i,j] = b*sqrt(ζ[j])*cos(θ[i])
                Y[i,j] = b*sqrt(ζ[j])*sin(θ[i])
                Z[i,j] = ζ[j]
                # X[i,j] = cos(θ[i])
                # Y[i,j] = sin(θ[i])
                # Z[i,j] = b*sqrt((ζ[j] - h1proj) / h)
            end
        end
    elseif typeof(s) == Hyperboloid
        R = s.R
        b = s.b
        for i in 1:nθ
            for j in 1:nζ
                X[i,j] = R*sqrt(1 + ζ[j]^2/b^2)*cos(θ[i])
                Y[i,j] = R*sqrt(1 + ζ[j]^2/b^2)*sin(θ[i])
                Z[i,j] = ζ[j]
            end
        end
    end

    Xc1 = [zeros(nθ, 1) X[:,1]]
    Xc2 = [zeros(nθ, 1) X[:,end]]
    Yc1 = [zeros(nθ, 1) Y[:,1]]
    Yc2 = [zeros(nθ, 1) Y[:,end]]
    Zc1 = [Z[:,1] Z[:,1]]
    Zc2 = [Z[:,1] Z[:,1]]

    # print(typeof(s))
    
    (X, Y, Z) = transform_to_axis(X, Y, Z, [0;0;1], a, ca1 - h1proj*a)
    (Xc1, Yc1, Zc1) = transform_to_axis(Xc1, Yc1, Zc1, [0;0;1], ps[1].a, ca1 + h1proj*a)   # FIX THESE: WRONG MESH POSITION
    (Xc2, Yc2, Zc2) = transform_to_axis(Xc2, Yc2, Zc2, [0;0;1], ps[2].a, ca2 + h1proj*a)

    return ((X,Xc1,Xc2), (Y,Yc1,Yc2), (Z,Zc1,Zc2))
end

function transform_to_axis(X, Y, Z, oldax, newax, center)
    # can do checks on sizes of X, Y, Z
    (nθ, nζ) = size(X)
    transform = zeros(3,3)
    if newax'*oldax == 1
        transform = I
    elseif newax'*oldax == -1
        transform = [1  0  0;
                     0  0 -1;
                     0 -1  0]
    else
        # cylinder axis oblique to z-axis:
        ax = cross(oldax, newax)
        ax = ax/norm(ax)
        cosang = newax'*oldax
        # Rodrigues's formula/axis-angle representation:
        axso3 = [ 0     -ax[3]   ax[2];
                  ax[3]  0      -ax[1];
                 -ax[2]  ax[1]   0]
        transform = I*cosang + (1 - cosang)*(ax*ax') + axso3*sqrt(1 - cosang^2)
    end

    # transform all cylinder points:
    for i = 1:nθ
        for j = 1:nζ
            newpos = transform*[X[i,j]; Y[i,j]; Z[i,j]]
            X[i,j] = newpos[1] #+ center[1]
            Y[i,j] = newpos[2] #+ center[2]
            Z[i,j] = newpos[3] #+ center[3]
        end
    end
    
    X .+= center[1]
    Y .+= center[2]
    Z .+= center[3]
    return (X, Y, Z)
end