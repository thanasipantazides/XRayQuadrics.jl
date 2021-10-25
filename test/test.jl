using XRayTrace
using LinearAlgebra

n = [0;0;1]
c = [1;1;1]

s = Plane(c,n)
qs = Quadric(s)

c = Cylinder(3, 4, [1;1;1], [1;0;0])
qc = Quadric(c)

p = Paraboloid(3, 4, [1;1;1], [0;0;1])
qp = Quadric(p)

h = Hyperboloid(3,4,8, [2;2;2], [0;1;0])
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