using XRayTrace
using LinearAlgebra
using GLMakie, Colors

n = [0;0;1]
c = [1;1;1]

s = Plane(c,n)
qs = Quadric(s)

c = Cylinder(3, [1;1;1], [1;0;0])
qc = Quadric(c)

p = Paraboloid(3, [1;1;1], [0;0;1])
qp = Quadric(p)

h = Hyperboloid(3,8, [2;2;2], [0;1;0])
qh = Quadric(h)

ray = Particle([0;0;10],[0;0;-1],10e3)

tp = inout(ray, qp)
tc = inout(ray, qc)
th = inout(ray, qh)

Ec = eigvals(qc.Q[1:3,1:3])
Eh = eigvals(qh.Q[1:3,1:3])
Ep = eigvals(qp.Q[1:3,1:3])
vc = eigvecs(qc.Q[1:3,1:3])
vh = eigvecs(qh.Q[1:3,1:3])
vp = eigvecs(qp.Q[1:3,1:3])

classify(qc)

ax = [1;1;1]
center = [12;10;0]
R = 1
b = 1
c = Cylinder(R, center, ax)
p = Paraboloid(R, center, ax)
h = Hyperboloid(R, b, center, ax)

p1 = Plane(c.a*1, ax)
p2 = Plane(c.a*4, ax)

tqc = TruncatedQuadric(c, [p1, p2])
tqp = TruncatedQuadric(p, [p1, p2])
tqh = TruncatedQuadric(h, [p1, p2])

f = Figure()
# scatter([0, 1],[0, 1],[0, 1], plot_color=:red)
# Makie.plot(tq)
# Makie.plot!(tq)

plot(tqp)
# plot!(tqh)