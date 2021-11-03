using LinearAlgebra

"""
    Plane(c::Vector{Float64}, a::Vector{Float64})

Construct a `Plane` from a point `c` which lies in the plane and a unit normal `a`. Currently, `c` and `a` are forced to be 3-vectors.
"""
struct Plane
    c::Vector{Float64}
    a::Vector{Float64}
    Plane(c,a) = begin
        if length(c) != 3 || length(a) != 3
            error("3-vectors please")
        else
            new(c,a/norm(a))
        end
    end
end

"""
    Cylinder(R::Float64, c::Vector{Float64}, a::Vector{Float64})

Construct a `Cylinder` with radius `R` and unit axis `a` which passes through point `c`.
"""
struct Cylinder
    R::Real          # radius
    c::Vector{Float64}  # center point
    a::Vector{Float64}  # axis (unit) (NOTE: c + h*a yields a point on the opposite cap => this is an inward normal axis)
    Cylinder(R,c,a) = begin
        if length(c) != 3 || length(a) != 3
            error("3-vectors please")
        else
            new(R,c,a/norm(a))
        end
    end
end

"""
    Paraboloid(R::Float64, c::Vector{Float64}, a::Vector{Float64})

Construct a `Paraboloid` at vertex `c` with growth rate `R` and unit axis `a`.
"""
struct Paraboloid
    R::Real
    c::Vector{Float64}
    a::Vector{Float64}
    Paraboloid(R,c,a) = begin
        if length(c) != 3 || length(a) != 3
            error("3-vectors please")
        else
            new(R,c,a/norm(a))
        end
    end
end

"""
    Hyperboloid(R::Float64, b::Float64, c::Vector{Float64}, a::Vector{Float64})

Construct a `Hyperboloid` of one sheet centered at `c` with minimum radius `R`, radial growth rate `b`, and unit axis `a`.
"""
struct Hyperboloid
    R::Real
    b::Real
    c::Vector{Float64}
    a::Vector{Float64}
    Hyperboloid(R,b,c,a) = begin
        if length(c) != 3 || length(a) != 3
            error("3-vectors please")
        else
            new(R,b,c,a/norm(a))
        end
    end
end


struct Quadric
    Q::Matrix{Float64}
    Quadric(Q) = begin
        if issymmetric(Q)
            new(Q)
        else
            error("quadric matrix must be symmetric")
        end
    end
end

"""
    TruncatedQuadric(q::Quadric, p::Vector{Plane})

Construct a `TruncatedQuadric` surface by slicing a non-degenerate `Quadric` q by two planes `p`. Planes may not contain quadric axis.
"""
struct TruncatedQuadric
    q::Quadric
    p::Vector{Plane}
    TruncatedQuadric(q,p) = begin
        Qh = q.Q[1:3,1:3]
        if all(eigvals(Qh) .!= 0)
            for plane in p
                if all(q.Q[1:3,1:3]*plane.a .== 0)
                    error("planes must be oblique to quadric axis")
                end
            end
        end
        new(q,p)
    end
end

function Quadric(s::Plane)
    Q = [zeros(3,3)      0.5*s.a;
         0.5*s.a'        -s.a'*s.c]
    return Quadric(Q)
end

function Quadric(s::Cylinder)
    Q = [s.a*s.a' - I            (I - s.a*s.a')*s.c;
         ((I - s.a*s.a')*s.c)'   s.R^2 + s.c'*(s.a*s.a' - I)*s.c]
    return Quadric(Q)
end

function Quadric(s::Paraboloid)
    Q = [s.a*s.a' - I   s.R^2/2*s.a + (I - s.a*s.a')*s.c;
         (s.R^2/2*s.a + (I - s.a*s.a')*s.c)'    -s.R^2*s.a'*s.c + s.c'*(s.a*s.a' - I)*s.c]
    return Quadric(Q)
end

function Quadric(s::Hyperboloid)
    γ = (1 + s.R^2/s.b^2)
    Q = [γ*s.a*s.a' - I     (I - γ*s.a*s.a')*s.c;
         ((I - γ*s.a*s.a')*s.c)'  s.R^2 + s.c'*(γ*s.a*s.a' - I)*s.c]
    return Quadric(Q)
end

TruncatedQuadric(s::Plane, p::Vector{Plane}) = TruncatedQuadric(Quadric(s), p) 
TruncatedQuadric(s::Cylinder, p::Vector{Plane}) = TruncatedQuadric(Quadric(s), p)
TruncatedQuadric(s::Paraboloid, p::Vector{Plane}) = TruncatedQuadric(Quadric(s), p)
TruncatedQuadric(s::Hyperboloid, p::Vector{Plane}) = TruncatedQuadric(Quadric(s), p)

"""
    TruncatedQuadric(q::Quadric, h1::Vector{Float64}, h2::Vector{Float64})

A convenience constructor for a `TruncatedQuadric`. Specify `h1` and `h2` as points along the axis of `q`, and the constructor will return `q` truncated by planes passing through `h1` and `h2` sharing a normal vector equal to the axis of `q`.
"""
function TruncatedQuadric(q::Quadric, h1::Vector{Float64}, h2::Vector{Float64})
    s = changerepresentation(q)
    if typeof(s) == Plane
        display("TruncatedQuadric just a bunch of planes")
    end
    
    if all((h1 - h2)/norm(h1 - h2) .== s.a/norm(s.a))
        p1 = Plane(h1, s.a)
        p2 = Plane(h2, s.a)
    else
        error("h1, h2 must lie on Quadric axis")
    end

    return TruncatedQuadric(q, [p1, p2])
end

function changerepresentation(q::Quadric)
    ε = 1e-15

    Qh = q.Q[1:3, 1:3]
    Qd = q.Q[1:3, end]
    Q0 = q.Q[end, end]
    E = eigvals(Qh)
    if all(E .== 0)
        # plane
        a = 2*Qd                    # Qh is zero, Qd is parallel to plane normal
        c = -Q0*a                   # Q0 = -a'*c: choose c to lie on a (free choice)
        return Plane(c, a)
        
    elseif all(abs.(E) .>= ε)
        # hyperboloid
        v = eigvecs(Qh)
        a = v[:,end]                # axis is eigenvector for largest eigenvalue
        c = (-Qh)\Qd                # center of hyperboloid
        γ = tr(Qh + I)              # Qh = γ*a*a' - I, and trace(a*a') == 1
        R = sqrt(Q0 - c'*Qh*c)      # Q0 = sum(c*c' .* Qh) + R^2
        b = 1/sqrt(γ - 1)*R         # γ = 1 + (R/b)^2
        return Hyperboloid(R, b, c, a)

    else
        v = eigvecs(Qh)
        a = v[:,end]        # axis is eigenvector for largest eigenvalue

        if rank(q.Q) == 4
            # paraboloid
            a = sign(Qd'*a)*a       # correct for antiparallel axis
            R = sqrt(2*Qd'*a)
            c = (-Qh)\(Qd - R^2/2*a)
            badI = findfirst(isnan.(c) .| isinf.(c))
            if !isnothing(badI)
                c[badI] = 1.0
            end
            return Paraboloid(R, c, a)
 
        elseif rank(q.Q) == 3
            # cylinder
            c = (-Qh)\Qd
            badI = findfirst(isnan.(c) .| isinf.(c))
            if !isnothing(badI)
                c[badI] = 1.0
            end
            R = sqrt(Q0 - c'*Qh*c)
            return Cylinder(R, c, a)
 
        else
            error("unclassified shape")

        end
    end
end

changerepresentation(s::Plane) = Quadric(s)
changerepresentation(s::Cylinder) = Quadric(s)
changerepresentation(s::Paraboloid) = Quadric(s)
changerepresentation(s::Hyperboloid) = Quadric(s)

function normal(q::Quadric, r::Vector{Float64})
    n = 2*q.Q*[r; 1]
    n = n[1:3]
    n = n/norm(n)
    return n
end

function normal(s::Plane, r::Vector{Float64})
    return s.a
end

function normal(s::Cylinder, r::Vector{Float64})
    n = 2*(s.a - (r - s.c)./norm(r - s.c))
    return n/norm(n)
end

function normal(s::Paraboloid, r::Vector{Float64})
    n = s.R^2*s.a + 2*(s.a - (r - s.c)./norm(r - s.c))
    return n/norm(n)
end

function normal(s::Hyperboloid, r::Vector{Float64})
    n = 2*((1 + s.R^2/s.b^2)*s.a - (r - s.c)./norm(r - s.c))
    return n/norm(n)
end

function classify(q::Quadric)

    Qh = q.Q[1:3,1:3]

    E = eigvals(Qh)

    p = sum(E .> 0)
    n = sum(E .< 0)
    z = sum(E .== 0)

    # options:
    #   [p, n, z]
    #   [0, 0, 3]:  plane
    #   [0, 1, 2]:  parabolic cylinder, parallel planes
    #   [0, 2, 1]:  imaginary cylinder
    #   [0, 3, 0]:  imag. ellipsoid
    #   [1, 0, 2]:  parabolic cylinder, parallel planes
    #   [1, 1, 1]:  hyperbolic paraboloid, hyperbolic cylinder
    #   [1, 2, 0]:  hyperboloid 2 sheet
    #   [2, 0, 1]:  elliptic paraboloid, elliptic cylinder
    #   [2, 1, 0]:  hyperboloid 1 sheet, cone @ origin
    #   [3, 0, 0]:  ellipsoid, imag. context
    shapes0 = [
        "none"  "none"  "none"  "ellipsoid";
        "none"  "none"  "hyperboloid"   "none";
        "none"  "hyperboloid"  "none"   "none";
        "ellipsoid" "none"  "none"  "none"
    ]
    shapes1 = [
        "none"  "none"  "elliptic cylinder/paraboloid"    "none";
        "none"  "hyperbolic paraboloid/cylinder"  "none"  "none";
        "elliptic cylinder/paraboloid" "none"   "none"  "none";
        "none"  "none"  "none"  "none"
    ]
    shapes2 = [
        "none"  "parabolic cylinder/planes" "none"  "none";
        "parabolic cylinder/planes" "none"  "none"  "none";
        "none"  "none"  "none"  "none";
        "none"  "none"  "none"  "none"
    ]
    shapes3 = [
        "plane"  "none"  "none"  "none";
        "none"   "none"  "none"  "none";
        "none"   "none"  "none"  "none";
        "none"   "none"  "none"  "none";
    ]

    shapes = cat(shapes0, shapes1, shapes2, shapes3, dims=3)

    shapename = shapes[p+1, n+1, z+1]
    return shapename
end