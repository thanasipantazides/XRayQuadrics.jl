using XRayTrace
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

@testset "XRayTrace.jl" begin
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
        @test tqp.R - p.R < ε
    end

    @testset "hyperboloid conversion" begin      # computed hyperbolic coeff and original coeff are same
        @test tqh.b - h.b < ε
    end

    @testset "intersections" begin
        @test all(in_out(sray, qs) .== (0,(s.a'*s.c - s.a'*sray.r0)/(s.a'*sray.v)))
        cray = Particle(c.c, [0;0;1], 10e3)
        @test all(in_out(cray, qc) .== (-c.R, c.R))
    end

    @testset "normal vectors" begin
        
    end
end
