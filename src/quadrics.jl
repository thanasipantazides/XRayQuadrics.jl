using LinearAlgebra

"""
    Plane(c::Vector{Float64}, a::Vector{Float64})

Construct a `Plane` from a point `c` which lies in the plane and a unit normal `a`.
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
    Paraboloid(b::Float64, c::Vector{Float64}, a::Vector{Float64})

Construct a `Paraboloid` at vertex `c` with growth rate `b` and unit axis `a`.
"""
struct Paraboloid
    b::Real
    c::Vector{Float64}
    a::Vector{Float64}
    Paraboloid(b,c,a) = begin
        if length(c) != 3 || length(a) != 3
            error("3-vectors please")
        else
            new(b,c,a/norm(a))
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

"""
    Quadric(Q::Matrix{Float64})

Construct a `Quadric` surface from a 4x4 symmetric quadric matrix.
"""
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
    TruncatedQuadric(q::Quadric, p::Vector{Plane}, caps)

Construct a `TruncatedQuadric` surface by slicing a non-reducible `Quadric` q with two planes `p`. Planes may not contain quadric axis. Argument `caps=true` closes the volume between the planes, `caps=false` defines an open-ended tube bounded by the `Quadric` and the planes.
"""
struct TruncatedQuadric
    q::Quadric
    p::Vector{Plane}
    caps::Bool
    TruncatedQuadric(q,p,caps) = begin
        Qr = q.Q[1:3,1:3]
        if all(eigvals(Qr) .!= 0)
            for plane in p
                if all(q.Q[1:3,1:3]*plane.a .== 0)
                    error("planes must be oblique to quadric axis")
                end
            end
        end
        new(q,p,caps)
    end
end

@doc raw"""
    Quadric(s::Plane)

Construct a `Quadric` from a `Plane` (with attributes axis = `a` = ``\boldsymbol{a}``; and point in plane = `c` = ``\boldsymbol{c}``) using the definition:

``\boldsymbol{Q}_s = \begin{bmatrix} \boldsymbol{0} & \frac{1}{2}\boldsymbol{a} \\ \frac{1}{2}\boldsymbol{a}^\mathsf{T} & -\boldsymbol{a}^\mathsf{T}\boldsymbol{c} \end{bmatrix}``
"""
function Quadric(s::Plane)
    Q = [zeros(3,3)      0.5*s.a;
         0.5*s.a'        -s.a'*s.c]
    return Quadric(Q)
end

@doc raw"""
    Quadric(s::Cylinder)

Construct a `Quadric` from a `Cylinder` (which has attributes axis = `a` = ``\boldsymbol{a}``; point on axis = `c` = ``\boldsymbol{c}``; and radius = `R` = ``R``) using the definitions:

``\boldsymbol{Q}_r = \boldsymbol{a}\boldsymbol{a}^\mathsf{T} - \boldsymbol{1}``

``\boldsymbol{Q}_c = \begin{bmatrix} \boldsymbol{Q}_r & -\boldsymbol{Q}_r\boldsymbol{c} \\ -\boldsymbol{c}^\mathsf{T}\boldsymbol{Q}_r & R^2 + \boldsymbol{c}^\mathsf{T}\boldsymbol{Q}_r\boldsymbol{c} \end{bmatrix}``
"""
function Quadric(s::Cylinder)
    Q = [s.a*s.a' - I            (I - s.a*s.a')*s.c;
         ((I - s.a*s.a')*s.c)'   s.R^2 + s.c'*(s.a*s.a' - I)*s.c]
    return Quadric(Q)
end

@doc raw"""
    Quadric(s::Paraboloid)

Construct a `Quadric` from a `Paraboloid` (which has attributes axis = `a` = ``\boldsymbol{a}``; point on axis = `c` = ``\boldsymbol{c}``; and quadratic parameter = `b` = ``b``) using the definitions:

``\boldsymbol{Q}_r = \boldsymbol{a}\boldsymbol{a}^\mathsf{T} - \boldsymbol{1}``

``\boldsymbol{Q}_p = \begin{bmatrix} \boldsymbol{Q}_r & \frac{b^2}{2}\boldsymbol{a} - \boldsymbol{Q}_r\boldsymbol{c} \\ \frac{b^2}{2}\boldsymbol{a}^\mathsf{T} - \boldsymbol{c}^\mathsf{T}\boldsymbol{Q}_r & -b^2\boldsymbol{a}^\mathsf{T}\boldsymbol{c} + \boldsymbol{c}^\mathsf{T}\boldsymbol{Q}_r\boldsymbol{c} \end{bmatrix}``
"""
function Quadric(s::Paraboloid)
    Q = [s.a*s.a' - I   s.b^2/2*s.a + (I - s.a*s.a')*s.c;
         (s.b^2/2*s.a + (I - s.a*s.a')*s.c)'    -s.b^2*s.a'*s.c + s.c'*(s.a*s.a' - I)*s.c]
    return Quadric(Q)
end

@doc raw"""
    Quadric(s::Hyperboloid)

Construct a `Quadric` from a `Hyperboloid` (which has attributes axis = `a` = ``\boldsymbol{a}``; point on axis = `c` = ``\boldsymbol{c}``; neck radius = `R` = ``R``; and slope parameter = `b` = ``b``) using the definitions:

``\gamma = 1 + \frac{R^2}{b^2}``

``\boldsymbol{Q}_r = \gamma\boldsymbol{a}\boldsymbol{a}^\mathsf{T} - \boldsymbol{1}``

``\boldsymbol{Q}_h = \begin{bmatrix} \boldsymbol{Q}_r & -\boldsymbol{Q}_r\boldsymbol{c} \\ -\boldsymbol{c}^\mathsf{T}\boldsymbol{Q}_r & R^2 + \boldsymbol{c}^\mathsf{T}\boldsymbol{Q}_r\boldsymbol{c} \end{bmatrix}``
"""
function Quadric(s::Hyperboloid)
    γ = (1 + s.R^2/s.b^2)
    Q = [γ*s.a*s.a' - I     (I - γ*s.a*s.a')*s.c;
         ((I - γ*s.a*s.a')*s.c)'  s.R^2 + s.c'*(γ*s.a*s.a' - I)*s.c]
    return Quadric(Q)
end

"""
    TruncatedQuadric(s::Cylinder, p::Vector{Plane}, caps=true)

Construct a `TruncatedQuadric` surface from a `Cylinder` and set of `Plane`s by converting `Cylinder` directly to a `Quadric`.
"""
TruncatedQuadric(s::Cylinder, p::Vector{Plane}, caps=true) = TruncatedQuadric(Quadric(s), p, caps)

"""
    TruncatedQuadric(s::Paraboloid, p::Vector{Plane}, caps=true)

Construct a `TruncatedQuadric` surface from a `Paraboloid` and set of `Plane`s by converting `Paraboloid` directly to a `Quadric`.
"""
TruncatedQuadric(s::Paraboloid, p::Vector{Plane}, caps=true) = TruncatedQuadric(Quadric(s), p, caps)

"""
    TruncatedQuadric(s::Hyperboloid, p::Vector{Plane}, caps=true)

Construct a `TruncatedQuadric` surface from a `Hyperboloid` and set of `Plane`s by converting `Hyperboloid` directly to a `Quadric`.
"""
TruncatedQuadric(s::Hyperboloid, p::Vector{Plane}, caps=true) = TruncatedQuadric(Quadric(s), p, caps)

"""
    TruncatedQuadric(q::Quadric, h1::Vector{Float64}, h2::Vector{Float64}, caps=true)

A convenience constructor for a `TruncatedQuadric`. Specify `h1` and `h2` as points along the axis of `q`, and the constructor will return `q` truncated by planes passing through `h1` and `h2` sharing a normal vector equal to the axis of `q`.
"""
function TruncatedQuadric(q::Quadric, h1::Vector{Float64}, h2::Vector{Float64}, caps=true)
    s = changerepresentation(q)
    if typeof(s) == Plane
        display("TruncatedQuadric just a bunch of planes")
    end
    
    if abs(((h1 - h2)/norm(h1 - h2))'*s.a/norm(s.a)) == 1
        p1 = Plane(h1, s.a)
        p2 = Plane(h2, s.a)
    else
        error("h1, h2 must lie on Quadric axis")
    end

    return TruncatedQuadric(q, [p1, p2], caps)
end

function shapepotential(s::Plane, r::Vector{Float64})
    return s.a'*r + s.c'*r
end

function shapepotential(s::Cylinder, r::Vector{Float64})
    return s.R^2 + (s.a'*(r - s.c))^2 - (r - s.c)'*(r - s.c)
end

function shapepotential(s::Paraboloid, r::Vector{Float64})
    return s.b^2*(s.a'*(r - s.c)) + (s.a'*(r - s.c))^2 - (r - s.c)'*(r - s.c)
end

function shapepotential(s::Hyperboloid, r::Vector{Float64})
    return s.R^2*(1 + 1/s.b^2*(s.a'*(r - s.c))^2) + (s.a'*(r - s.c))^2 - (r - s.c)'*(r - s.c)
end

function shapepotential(q::Quadric, r::Vector{Float64})
    return [r; 1]'*q.Q*[r; 1]
end

function normal_grad(q::Union{Quadric, Plane, Cylinder, Hyperboloid, Paraboloid}, r::Vector{Float64})
    n = gradient(p -> shapepotential(q, p), r)[1]
    return n/norm(n)
end

@doc raw"""
    changerepresentation(q::Quadric)

Converts a `Quadric` to a `Plane`, `Cylinder`, `Paraboloid`, or `Hyperboloid` based on its matrix structure. 

!!! note 
    This function is not fast. Minimize conversion between explicit surface and quadric matrix in heavy-lifting code. 
"""
function changerepresentation(q::Quadric)
    ε = 1e-15

    Qr = q.Q[1:3, 1:3]
    qd = q.Q[1:3, end]
    q0 = q.Q[end, end]
    E = eigvals(Qr)
    if all(E .== 0)
        # plane
        a = 2*qd                    # Qr is zero, qd is parallel to plane normal
        c = -q0*a                   # q0 = -a'*c: choose c to lie on a (free choice)
        return Plane(c, a)
        
    elseif all(abs.(E) .>= ε)
        # hyperboloid
        v = eigvecs(Qr)
        a = v[:,end]                # axis is eigenvector for largest eigenvalue
        c = (-Qr)\qd                # center of hyperboloid
        γ = tr(Qr + I)              # Qr = γ*a*a' - I, and trace(a*a') == 1
        R = sqrt(q0 - c'*Qr*c)      # q0 = sum(c*c' .* Qr) + R^2
        b = 1/sqrt(γ - 1)*R         # γ = 1 + (R/b)^2
        return Hyperboloid(R, b, c, a)

    else
        v = eigvecs(Qr)
        a = v[:,end]        # axis is eigenvector for largest eigenvalue

        if rank(q.Q) == 4
            # paraboloid
            a = sign(qd'*a)*a       # correct for antiparallel axis
            b = sqrt(2*qd'*a)
            
            cb1 = v[:,1]'*qd
            cb2 = v[:,2]'*qd
            cb3 = -(q0 + cb1^2 + cb2^2)/b^2

            c = v*[cb1; cb2; cb3]

            return Paraboloid(b, c, a)
 
        elseif rank(q.Q) == 3
            # cylinder
            c = (-Qr)\qd
            badI = findfirst(isnan.(c) .| isinf.(c))
            if !isnothing(badI)
                c[badI] = 1.0
            end
            R = sqrt(q0 - c'*Qr*c)
            return Cylinder(R, c, a)
 
        else
            error("unclassified shape")

        end
    end
end

"""
    changerepresentation(s::Plane)

Convert a `Plane` to `Quadric` representation via the `Quadric(s::Plane)` constructor.
"""
changerepresentation(s::Plane) = Quadric(s)

"""
    changerepresentation(s::Cylinder)

Convert a `Cylinder` to `Quadric` representation via the `Quadric(s::Cylinder)` constructor.
"""
changerepresentation(s::Cylinder) = Quadric(s)

"""
    changerepresentation(s::Paraboloid)

Convert a `Paraboloid` to `Quadric` representation via the `Quadric(s::Paraboloid)` constructor.
"""
changerepresentation(s::Paraboloid) = Quadric(s)

"""
    changerepresentation(s::Hyperboloid)

Convert a `Hyperboloid` to `Quadric` representation via the `Quadric(s::Hyperboloid)` constructor.
"""
changerepresentation(s::Hyperboloid) = Quadric(s)

"""
    normal(q::Quadric, r::Vector{Float64})

Returns unit normal vector to a `Quadric` at position `r` on the surface. 

!!! warning
    This will return a normal vector even if `r` does not fall on the surface. It is up to *you* to input an `r` you know to lie on the surface of `q`.
"""
function normal(q::Quadric, r::Vector{Float64})
    n = 2*q.Q*[r; 1]
    n = n[1:3]
    n = n/norm(n)
    return n
end

"""
    normal(s::Plane, r::Vector{Float64})

Returns unit normal vector to a `Plane`. This will return the `Plane`'s normal vector regardless of the `r` supplied; it takes the `r` argument nonetheless to maintain the same pattern as `normal()` functions for nonreducible `Quadric`s.
"""
function normal(s::Plane, r::Vector{Float64})
    return s.a
end

"""
    normal(s::Cylinder, r::Vector{Float64})

Returns unit normal vector to a `Cylinder` at position `r` on the surface. 

!!! warning
    This will return a normal vector even if `r` does not fall on the surface. It is up to *you* to input an `r` you know to lie on the surface of `s`.
"""
function normal(s::Cylinder, r::Vector{Float64})
    n = 2*(s.a - (r - s.c)./norm(r - s.c))
    return n/norm(n)
end

"""
    normal(s::Paraboloid, r::Vector{Float64})

Returns unit normal vector to a `Paraboloid` at position `r` on the surface. 

!!! warning
    This will return a normal vector even if `r` does not fall on the surface. It is up to *you* to input an `r` you know to lie on the surface of `s`.
"""
function normal(s::Paraboloid, r::Vector{Float64})
    n = s.b^2*s.a + 2*(s.a - (r - s.c)./norm(r - s.c))
    return n/norm(n)
end

"""
    normal(s::Hyperboloid, r::Vector{Float64})

Returns normal vector to a `Hyperboloid` at position `r` on the surface. 

!!! warning
    This will return a normal vector even if `r` does not fall on the surface. It is up to *you* to input an `r` you know to lie on the surface of `q`.
"""
function normal(s::Hyperboloid, r::Vector{Float64})
    n = 2*((1 + s.R^2/s.b^2)*s.a - (r - s.c)./norm(r - s.c))
    return n/norm(n)
end

function classify(q::Quadric)

    Qr = q.Q[1:3,1:3]

    E = eigvals(Qr)

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