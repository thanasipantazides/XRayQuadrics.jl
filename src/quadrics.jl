using LinearAlgebra

struct Plane
    c::Vector{Float64}
    a::Vector{Float64}
    Plane(c,a) = begin
        if length(c) != 3 || length(a) != 3
            error("out of order")
        else
            new(c,a/norm(a))
        end
    end
end

struct Cylinder
    R::Real          # radius
    c::Vector{Float64}  # center point
    a::Vector{Float64}  # axis (unit) (NOTE: c + h*a yields a point on the opposite cap => this is an inward normal axis)
    Cylinder(R,c,a) = begin
        if length(c) != 3 || length(a) != 3
            error("out of order")
        else
            new(R,c,a/norm(a))
        end
    end
end

# NEED SOME WAY TO DISPLACE THESE FROM c: h1 and h2?
struct Paraboloid
    R::Real
    c::Vector{Float64}
    a::Vector{Float64}
    Paraboloid(R,c,a) = begin
        if length(c) != 3 || length(a) != 3
            error("out of order")
        else
            new(R,c,a/norm(a))
        end
    end
end

struct Hyperboloid
    R::Real
    b::Real
    c::Vector{Float64}
    a::Vector{Float64}
    Hyperboloid(R,b,c,a) = begin
        if length(c) != 3 || length(a) != 3
            error("out of order")
        else
            new(R,b,c,a/norm(a))
        end
    end
end

struct Quadric
    Q::Matrix{Float64}
end

struct TruncatedQuadric
    Q::Quadric
    P::Vector{Plane}
end

# quadric constructors
function Quadric(Q::Matrix{Real})
    # do some definiteness checks or something
    return Quadric(Q)
end

# function TruncatedQuadric(Q::Quadric, P::Vector{Plane})

# end

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
         (s.R^2/2*s.a + (I - s.a*s.a')*s.c)'    s.R^2*s.a'*s.c + s.c'*(s.a*s.a' - I)*s.c]
    return Quadric(Q)
end

function Quadric(s::Hyperboloid)
    γ = (1 + s.R^2/s.b^2)
    Q = [γ*s.a*s.a' - I     (I - γ*s.a*s.a')*s.c;
         ((I - γ*s.a*s.a')*s.c)'  s.R^2 + s.c'*(γ*s.a*s.a' - I)*s.c]
    return Quadric(Q)
end

function changerepresentation(q::Quadric)
    Qh = q.Q[1:3, 1:3]
    Qd = q.Q[1:3, end]
    Q0 = q.Q[end, end]
    E = eigvals(Qh)
    if all(E .== 0)
        # plane
        a = 2*Qd                        # Qh is zero, Qd is parallel to plane normal
        # a = a/norm(a)                   
        c = -Q0*a                       # Q0 = -a'*c: choose c to lie on a (free choice)
        return Plane(c, a)
        
    elseif all(E .!= 0)
        # hyperboloid
        v = eigvecs(Qh)
        a = v[:,end]                    # axis is eigenvector for largest eigenvalue
        # a = a/norm(a)
        c = (-Qh)\Qd                 # center of hyperboloid
        γ = tr(Qh + I)               # Qh = γ*a*a' - I, and trace(a*a') == 1
        R = sqrt(Q0 - c'*Qh*c)  # Q0 = sum(c*c' .* Qh) + R^2
        b = 1/sqrt(γ - 1)*R             # γ = 1 + (R/b)^2
        return Hyperboloid(R, b, c, a)

    else
        v = eigvecs(Qh)
        a = v[:,end]        # axis is eigenvector for largest eigenvalue
        # a = a/norm(a)

        if rank(q.Q) == 4
            # paraboloid
            a = sign(Qd'*a)*a       # correct for antiparallel axis
            R = sqrt(2*Qd'*a)
            c = (-Qh)\(Qd - R^2/2*a)
            badI = findfirst(isnan.(c) .| isinf.(c))
            c[badI] = 1.0
            return Paraboloid(R, c, a)
        elseif rank(q.Q) == 3
            # cylinder
            c = (-Qh)\Qd
            badI = findfirst(isnan.(c) .| isinf.(c))
            c[badI] = 1.0
            R = sqrt(Q0 - c'*Qh*c)
            return Cylinder(R, c, a)
        else
            error("unclassified shape")
        end
    end
end

function changerepresentation(s::Plane)
    return Quadric(s)
end

function changerepresentation(s::Cylinder)
    return Quadric(s)
end

function changerepresentation(s::Paraboloid)
    return Quadric(s)
end

function changerepresentation(s::Hyperboloid)
    return Quadric(s)
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