using XRayQuadrics
using LinearAlgebra
using Colors
using Plots

axis = [0;0;1]
center = [0;0;0]

R = 1
b = 1

s = Plane(center, axis)
qs = Quadric(s)

c = Cylinder(R, center, axis)
qc = Quadric(c)

p = Paraboloid(R, center, axis)
qp = Quadric(p)

h = Hyperboloid(R,b, center, axis)
qh = Quadric(h)

p1 = Plane(c.a*0.01, axis)
p2 = Plane(c.a*2, axis)

tqc = TruncatedQuadric(c, [p1, p2], true)
tqp = TruncatedQuadric(p, [p1, p2], true)
tqh = TruncatedQuadric(h, [p1, p2], true)

gr()
splot = Plots.surface(tqp, color=:blue, plot_caps=true, size=(500,500))
display(splot)
Plots.surface!(tqh, color=:red)
Plots.surface!(tqc, color=:green)
Plots.plot!([center center+3*axis][1,:], [center center+3*axis][2,:], [center center+3*axis][3,:], color=:blue, linewidth=4)