using XRayQuadrics
using LinearAlgebra
using Test

s = Plane([0;0;1],[1;1;1])
qs = Quadric(s)
tqs = changerepresentation(qs)

c = Cylinder(3, [1;1;1], [1;0;0])
qc = Quadric(c)
tqc = changerepresentation(qc)

h = Hyperboloid(3, 8, [2;2;2], [0;1;0])
qh = Quadric(h)
tqh = changerepresentation(qh)

p = Paraboloid(3, [3;3;3], [0;0;1])
qp = Quadric(p)
tqp = changerepresentation(qp)

sray = Particle([0;0;0],[1;1;1],10e3)

n = 10
solA = rand(n,n)
solb = rand(n)
solc = 10

@testset "XRayQuadrics.jl" begin
    # Write your tests here.
    # test boolean expressions with @test
    ε = 1e-15
    @testset "axis conversion" begin
        @test norm(tqs.a'*s.a) - 1 < ε      # computed axis and original axis are parallel or antiparallel
        @test norm(tqc.a'*c.a) - 1 < ε
        @test norm(tqh.a'*h.a) - 1 < ε
        @test norm(tqp.a'*p.a) - 1 < ε
    end

    @testset "center conversion" begin
        @test s.a'*(tqs.c - s.c) < ε        # computed center and original center are coplanar
        @test c.a'*(tqc.c - c.c) - 1 < ε    # computed and true centers lie on axis
        @test norm(tqp.c - p.c) < ε         # computed and true centers are same
        @test norm(tqh.c - h.c) < ε
    end

    @testset "radius conversion" begin
        @test tqc.R - c.R < ε               # computed radius and original radius are same
        @test tqh.R - h.R < ε
    end

    @testset "hyperboloid, paraboloid conversion" begin      # computed hyperbolic coeff and original coeff are same
        @test tqh.b - h.b < ε
        @test tqp.b - p.b < ε
    end

    @testset "intersections" begin
        @test all(in_out(sray, qs) .== ((s.a'*s.c - s.a'*sray.r0)/(s.a'*sray.v),0))
        cray = Particle(c.c, [0;0;1], 10e3)
        @test all(in_out(cray, qc) .== (-c.R, c.R))
    end

    @testset "normal vectors" begin
        
    end

    @testset "quadratic solution" begin
        cyl = Cylinder(1, [1;1;1], [0;0;1])
        qcyl = Quadric(cyl)
        
        cylA = qcyl.Q[1:3,1:3]
        cylb = qcyl.Q[1:3,4]
        cylc = qcyl.Q[4,4]

        r0 = [-10;1;1]
        v0 = [1;0;0]

        V = [v0 zeros(3,1); 
                zeros(3,1) v0]
        p = [r0;r0]
        Q = [cylA zeros(3,3);
        zeros(3,3) cylA]
        solA = V'*Q*V
        solb = 2*Q*V

        solc = p'*Q*p - 2*cylc
        (xint1, xint2) = solve_quadratic(solA,solb,solc)

        println(xint1, xint2)
        
        # (x1, x2) = solve_quadratic(solA,solb,solc)
    end
end
