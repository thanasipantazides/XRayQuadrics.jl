struct Quadric
    Q::Matrix{Float64}
end

struct Plane
    c::Vector{Float64}
    n::Vector{Float64}
end

struct Cylinder
    R::Float64          # radius
    h::Float64          # height
    c::Vector{Float64}  # center point
    a::Vector{Float64}  # axis (unit) (NOTE: c + h*a yields a point on the opposite cap => this is an inward normal axis)
end

# NEED SOME WAY TO DISPLACE THESE FROM c: h1 and h2?
struct Paraboloid
    R::Float64
    h::Float64
    c::Vector{Float64}
    a::Vector{Float64}
end

struct Hyperboloid
    R::Float64
    h::Float64
    b::Float64
    c::Vector{Float64}
    a::Vector{Float64}
end

function Plane(c::Vector{Real}, n::Vector{Real})
    if length(n) != 3 || length(c) != 3
        error("3-vectors, please")
    end

    return Plane(c, n/norm(n))
end

function Cylinder(R::Real, h::Real, c::Vector{Real}, a::Vector{Real})
    if length(a) != 3 || length(c) != 3
        error("3-vectors, please")
    end

    return Cylinder(R, h, c, a/norm(a))
end

function Paraboloid(R::Real, h::Real, c::Vector{Real}, a::Vector{Real})
    if length(a) != 3 || length(c) != 3
        error("3-vectors, please")
    end

    return Paraboloid(R, h, c, a/norm(a))
end

function Hyperboloid(R::Real, h::Real, b::Real, c::Vector{Real}, a::Vector{Real})
    if length(a) != 3 || length(c) != 3
        error("3-vectors, please")
    end

    return Hyperboloid(R, h, b, c, a/norm(a))
end

# quadric constructors
function Quadric(Q::Matrix{Real})
    # do some definiteness checks or something
    return Quadric(Q)
end

function Quadric(s::Plane)
    Q = [zeros(3,3)      0.5*s.n;
         0.5*s.n'        -s.n'*s.c]
    return Quadric(Q)
end

function Quadric(s::Cylinder)
    Q = [s.a*s.a' - I            (I - s.a*s.a')*s.c;
         ((I - s.a*s.a')*s.c)'   s.R^2 + sum(s.c*s.c' .* (s.a*s.a' - I))]
    return Quadric(Q)
end

function Quadric(s::Paraboloid)
    Q = [s.a*s.a' - I   s.R^2/2*s.a + (I - s.a*s.a')*s.c;
         (s.R^2/2*s.a + (I - s.a*s.a')*s.c)'    s.R^2*s.a'*s.c + sum(s.c*s.c' .* (s.a*s.a' - I))]
    return Quadric(Q)
end

function Quadric(s::Hyperboloid)
    g = (1 + s.R^2/s.b^2)
    Q = [g*s.a*s.a' - I     (I - s.a*s.a'*g)*s.c;
         ((I - s.a*s.a'*g)*c)'  s.R^2 + sum(s.c*s.c' .* (s.a*s.a'*g - I))]
    return Quadric(Q)
end

function inout(p::Particle, q::Quadric)
    # Project to homogeneous coordinate
    rh = [p.r0;1]
    vh = [p.v;1]
    Q = q.Q
    # check for coplanar/axial alignment
    denom = vh'*Q*vh
    if denom == 0.0
        return (0,0)
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