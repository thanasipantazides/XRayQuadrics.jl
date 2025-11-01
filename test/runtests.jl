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

sray = Particle([0;0;0],[1;1;1],10e3, 0)

n = 10
solA = rand(n,n)
solb = rand(n)
solc = 10

@testset "XRayQuadrics.jl" begin
    # Write your tests here.
    # test boolean expressions with @test
    ε = 1e-12 # note: changerepresentation tolerance is 1e-12
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

    # @testset "hyperboloid conversion" begin      # computed hyperbolic coeff and original coeff are same
    #     @test tqh.b - h.b < ε
    # end
    
    @testset "cone conversion" begin
        np = 100
        inh = zeros(np)
        outh = zeros(np)
        ina = zeros(3, np)
        outa = zeros(3, np)
        inc = zeros(3, np)
        outc = zeros(3, np)
        for k = 1:np
            ax = rand(3)*2 .- 1
            ax = ax/norm(ax)
            cent = rand(3)*4 .- 2
            cone = Cone(
                rand()*3 + 0.5, 
                cent, 
                ax
            )
            qcone = Quadric(cone)
            recone = changerepresentation(qcone)
            
            inh[k] = cone.h
            outh[k] = recone.h
            ina[:,k] = cone.a
            outa[:,k] = recone.a
            inc[:,k] = cone.c
            outc[:,k] = recone.c
        end
        @testset "h" begin [@test abs(inh[k] - outh[k]) < ε for k in 1:np] end
        @testset "a" begin [@test abs(ina[:,k]'*outa[:,k] - 1) < ε for k in 1:np] end
        @testset "|a'*a| = 1" begin [@test abs(ina[:,k]'*outa[:,k]) - 1 < ε for k in 1:np] end
        @testset "c" begin [@test abs(sum(inc[:,k] .- outc[:,k])) < 3*ε for k in 1:np] end
    end
    
    @testset "cylinder conversion" begin
        np = 100
        inr = zeros(np)
        outr = zeros(np)
        ina = zeros(3, np)
        outa = zeros(3, np)
        inc = zeros(3, np)
        outc = zeros(3, np)
        for k = 1:np
            ax = rand(3)*2 .- 1
            ax = ax/norm(ax)
            cent = rand(3)*4 .- 2
            cyl = Cylinder(
                rand()*2 + 0.1, 
                cent, 
                ax
            )
            qcyl = Quadric(cyl)
            recyl = changerepresentation(qcyl)
            
            inr[k] = cyl.R
            outr[k] = recyl.R
            ina[:,k] = cyl.a
            outa[:,k] = recyl.a
            inc[:,k] = cyl.c
            outc[:,k] = recyl.c
        end
        @testset "R" begin [@test abs(inr[k] - outr[k]) < ε for k in 1:np] end
        @testset "a" begin [@test abs(ina[:,k]'*outa[:,k] - 1) < ε for k in 1:np] end
        @testset "|a'*a| = 1" begin [@test abs(ina[:,k]'*outa[:,k]) - 1 < ε for k in 1:np] end
        @testset "c" begin [@test abs(sum(inc[:,k] .- outc[:,k])) < 3*ε for k in 1:np] end
    end
    
    @testset "paraboloid conversion" begin
        np = 100
        inb = zeros(np)
        outb = zeros(np)
        ina = zeros(3, np)
        outa = zeros(3, np)
        inc = zeros(3, np)
        outc = zeros(3, np)
        for k = 1:np
            ax = rand(3)*2 .- 1
            ax = ax/norm(ax)
            cent = rand(3)*4 .- 2
            par = Paraboloid(
                rand()*5, 
                cent, 
                ax
            )
            qpar = Quadric(par)
            repar = changerepresentation(qpar)
            
            inb[k] = par.b
            outb[k] = repar.b
            ina[:,k] = par.a
            outa[:,k] = repar.a
            inc[:,k] = par.c
            outc[:,k] = repar.c
        end
        @testset "b" begin [@test abs(inb[k] - outb[k]) < ε for k in 1:np] end
        @testset "a" begin [@test abs(ina[:,k]'*outa[:,k] - 1) < ε for k in 1:np] end
        @testset "|a'*a| = 1" begin [@test abs(ina[:,k]'*outa[:,k]) - 1 < ε for k in 1:np] end
        @testset "c" begin [@test abs(sum(inc[:,k] .- outc[:,k])) < 3*ε for k in 1:np] end
    end
    
    @testset "ellipsoid conversion" begin
        np = 100
        ine = zeros(np)
        oute = zeros(np)
        ind = zeros(np)
        outd = zeros(np)
        ina = zeros(3, np)
        outa = zeros(3, np)
        inc = zeros(3, np)
        outc = zeros(3, np)
        for k = 1:np
            ax = rand(3)*2 .- 1
            ax = ax/norm(ax)
            cent = rand(3)*4 .- 2
            d = rand()*4 - 1
            e = d - rand()
            ell = Ellipsoid(
                e,
                d,
                cent, 
                ax
            )
            qell = Quadric(ell)
            reell = changerepresentation(qell)
            
            ind[k] = ell.d
            outd[k] = reell.d
            ine[k] = ell.e
            oute[k] = reell.e
            ina[:,k] = ell.a
            outa[:,k] = reell.a
            inc[:,k] = ell.c
            outc[:,k] = reell.c
        end
        @testset "d (semimajor)" begin [@test abs(ind[k] - outd[k]) < ε for k in 1:np] end
        @testset "e (semiminor)" begin [@test abs(ine[k] - oute[k]) < ε for k in 1:np] end
        @testset "a" begin [@test abs(ina[:,k]'*outa[:,k] - 1) < ε for k in 1:np] end
        @testset "|a'*a| = 1" begin [@test abs(ina[:,k]'*outa[:,k]) - 1 < ε for k in 1:np] end
        @testset "c" begin [@test abs(sum(inc[:,k] .- outc[:,k])) < 3*ε for k in 1:np] end
    end

    @testset "hyperboloid conversion" begin
        np = 100
        inr = zeros(np)
        outr = zeros(np)
        inb = zeros(np)
        outb = zeros(np)
        ina = zeros(3, np)
        outa = zeros(3, np)
        inc = zeros(3, np)
        outc = zeros(3, np)
        for k = 1:np
            ax = rand(3)*2 .- 1
            ax = ax/norm(ax)
            cent = rand(3)*4 .- 2
            R = rand()*4 - 1
            b = rand()*3 - 1
            hyp = Hyperboloid(
                R,
                b,
                cent, 
                ax
            )
            qhyp = Quadric(hyp)
            rehyp = changerepresentation(qhyp)
            
            inr[k] = hyp.R
            outr[k] = rehyp.R
            inb[k] = hyp.b
            outb[k] = rehyp.b
            ina[:,k] = hyp.a
            outa[:,k] = rehyp.a
            inc[:,k] = hyp.c
            outc[:,k] = rehyp.c
        end
        @testset "R" begin [@test abs(inr[k] - outr[k]) < ε for k in 1:np] end
        @testset "b" begin [@test abs(inb[k] - outb[k]) < ε for k in 1:np] end
        @testset "a" begin [@test abs(ina[:,k]'*outa[:,k] - 1) < ε for k in 1:np] end
        @testset "|a'*a| = 1" begin [@test abs(ina[:,k]'*outa[:,k]) - 1 < ε for k in 1:np] end
        @testset "c" begin [@test abs(sum(inc[:,k] .- outc[:,k])) < 3*ε for k in 1:np] end
    end
    
    @testset "intersections" begin
        @test all(in_out(sray, qs) .== ((s.a'*s.c - s.a'*sray.r0)/(s.a'*sray.v),0))
        cray = Particle(c.c, [0;0;1], 10e3, 1)
        @test all(in_out(cray, qc) .== (-c.R, c.R))
        
        @testset "parallel ray intersection" begin
            h = 1
            truncs = [
                Plane(c.c - c.a*h/2, c.a),
                Plane(c.c + c.a*h/2, c.a),
            ]
            tq = TruncatedQuadric(c, truncs, false)
            pdir = cross(c.a, rand(3))
            pdir = pdir/norm(pdir) # direction to move test points (perpendicular to axis)
            pos0 = c.c + (c.a'*c.c)*c.a/8 # initial position (along the axis, a little away from c)
            n = 100
            for k in 1:n
                p = pos0 + 2*c.R/n*(k - 1)*pdir
                res = inside(tq, p)
                if k <= n/2
                    @test res
                    
                else
                    @test !res
                end
            end
            
            ppar = Particle(
                pos0 - c.R*pdir, # displace along plane-parallel directoin
                pdir,            # fly parallel to planes
                1,
                2
            )
            pper = Particle(
                pos0 - c.R*c.a,  # displace along quadric axis
                c.a,            # fly parallel to quadric axis
                1,
                2
            )
            tpar = interactiontimes(tq, ppar)
            tper = interactiontimes(tq, pper)
            @test (tpar[2] - tpar[1]) - 2*c.R < ε
            @test (tper[2] - tper[1]) - h < ε
        end
        @testset "axis-plane" begin
            h = 1
            truncs = [
                Plane(c.c - c.a*h/2, c.a),
                Plane(c.c + c.a*h/2, c.a),
                # Plane(c.c - rand(3)*h*4, c.a),
                # Plane(c.c + rand(3)*h*4, c.a),
            ]
            tq = TruncatedQuadric(c, truncs, false)
            
            ca1 = axis_plane_intersection(c, truncs[1])
            ca2 = axis_plane_intersection(c, truncs[2])
            tqca1 = axis_plane_intersection(c, truncs[1])
            tqca2 = axis_plane_intersection(c, truncs[2])
            
            
            @test abs(sum(ca1 .- truncs[1].c)) < 3*ε
            @test abs(sum(ca2 .- truncs[2].c)) < 3*ε
            
            @test abs(sum(tqca1 .- truncs[1].c)) < 3*ε
            @test abs(sum(tqca2 .- truncs[2].c)) < 3*ε
        end
    end

    @testset "normal vectors" begin
        h = 1
        truncs = [
            Plane(c.c - c.a*h/2, c.a),
            Plane(c.c + c.a*h/2, c.a),
        ]
    end

    @testset "block properties" begin
        cyls = [
            Cylinder(1.0, [1;1;1], [0;0;1]),
            Cylinder(2.0, [1;1;1], [0;0;1]),
            Cylinder(3.0, [1;1;1], [0;0;1]),
            Cylinder(4.0, [1;1;1], [0;0;1]),
            Cylinder(5.0, [1;1;1], [0;0;1])
        ]
        n = length(cyls)
        qs = Vector{Quadric}(undef, n)
        
        bigQ = zeros(n*4, n*4)
        bigV = zeros(n*4, n)
        bigP = zeros(n*4, 1)
        v = [1;2;3;0]
        v = v/norm(v)
        p = [8;4;32;1]

        for (i,cyl) in enumerate(cyls)
            qs[i] = Quadric(cyl)
            bigI = (i - 1)*4 + 1
            bigQ[bigI:bigI+3, bigI:bigI+3] = qs[i].Q
            bigV[bigI:bigI+3, i] = v
            bigP[bigI:bigI+3] = p
        end

        quadblock = bigV'*bigQ*bigV
        linblock = 2*bigP'*bigQ*bigV
        
        # display(bigQ)
        
        # display(eigvals(bigQ))
        # println("Q rank deficiency: ",rank(bigQ) - 5*4)
        # println("Q pinv deficiency: ",rank(pinv(bigQ)) - 5*4)
        # println("V'QV")
        # display(quadblock)
        # display(pinv(quadblock))
        # println("V'QV rank deficiency: ",rank(quadblock) - 5)
        # println("V'QV pinv deficiency: ",rank(pinv(quadblock)) - 5)
        # println("(V'QV)^-1 QV")
        # display(pinv(quadblock)*linblock')

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
        # (xint1, xint2) = solve_quadratic(solA,solb,solc)

        # println(xint1, xint2)
        
        # (x1, x2) = solve_quadratic(solA,solb,solc)
    end
end
