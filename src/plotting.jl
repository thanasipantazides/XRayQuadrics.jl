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

function plot!(ax::GLMakie.LScene, tq::TruncatedQuadric; caps=true)
    (X, Y, Z, T) = cartesian_grid(tq)

    colors = Dict(Paraboloid=> :green, Hyperboloid=> :blue, Cylinder=> :red, Cone=> :yellow, Ellipsoid=> :magenta)

    surface!(
        ax, 
        X[1], Y[1], Z[1],
        colorrange=(-30, -20),
        highclip=(colors[T], 0.5),
        transparency=true
    )
    
    if caps
        surface!(
            ax, 
            X[2], Y[2], Z[2],
            colorrange=(-30, -20),
            highclip=(colors[T], 0.3),
            transparency=true
        )
        surface!(
            ax, 
            X[3], Y[3], Z[3],
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
    (Xs, Ys, Zs) = get_mesh(s, ps, nθ, nζ)
    return (Xs, Ys, Zs, T)
end

function get_mesh(s::Union{Cylinder, Paraboloid, Hyperboloid, Cone, Ellipsoid}, ps::Vector{Plane}, nθ::Int64, nζ::Int64)

    ca1 = axis_plane_intersection(s, ps[1])
    ca2 = axis_plane_intersection(s, ps[2])
    
    projaca1 = (ca1 - s.c)'*s.a
    projaca2 = (ca2 - s.c)'*s.a
    if projaca2 < projaca1
        temp = ca1
        ca1 = ca2
        ca2 = temp
    end

    h = norm(ca2 - ca1)

    h1 = ca1'*s.a
    h2 = ca2'*s.a

    a = s.a
    c = s.c
    θ = range(0, stop=2π, length=nθ)
    ζ = range(0, stop=h, length=nζ)
    # ζ = range(h1, stop=h2, length=nζ)

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
        h = s.h
        for i in 1:nθ
            for j in 1:nζ
                X[i,j] = 1/h*ζ[j]*cos(θ[i])
                Y[i,j] = 1/h*ζ[j]*sin(θ[i])
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
                X[i,j] = d*cos(θ[i])*sin(acos(ζ[j]/h))
                Y[i,j] = d*sin(θ[i])*sin(acos(ζ[j]/h))
                Z[i,j] = ζ[j]
            end
        end
    elseif typeof(s) == Paraboloid
        b = s.b
        for i in 1:nθ
            for j in 1:nζ
                X[i,j] = b*sqrt(ζ[j])*cos(θ[i])
                Y[i,j] = b*sqrt(ζ[j])*sin(θ[i])
                Z[i,j] = ζ[j]
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
    
    (X, Y, Z) = transform_to_axis(X, Y, Z, [0;0;1], a, ca1)
    (Xc1, Yc1, Zc1) = transform_to_axis(Xc1, Yc1, Zc1, [0;0;1], ps[1].a, ca1)   # FIX THESE: WRONG MESH POSITION
    (Xc2, Yc2, Zc2) = transform_to_axis(Xc2, Yc2, Zc2, [0;0;1], ps[2].a, ca2)

    return ((X,Xc1,Xc2), (Y,Yc1,Yc2), (Z,Zc1,Zc2))
end

function axis_plane_intersection(s::Union{Cylinder, Paraboloid, Hyperboloid, Cone, Ellipsoid}, p::Plane)
    # project plane points onto axis
    # ca = s.c + (p.a'*p.c - p.a'*s.c)/(p.a'*s.a)*s.a
    ca = s.a'*(p.c - s.c)*s.a + s.c
    return ca
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