using XRayTrace
using LinearAlgebra
using GLMakie, Colors
using Plots

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

ax = [1;2;3]
ax = ax/norm(ax)
center = [10;0;0]
R = 1
b = 1
c = Cylinder(R, center, ax)
p = Paraboloid(R, center, ax)
h = Hyperboloid(R, b, center, ax)

p1 = Plane(c.a*0.01, ax)
p2 = Plane(c.a*6, ax)

qc = Quadric(c)
qp = Quadric(p)
qh = Quadric(h)

tqc = TruncatedQuadric(c, [p1, p2])
tqp = TruncatedQuadric(p, [p1, p2])
tqh = TruncatedQuadric(h, [p1, p2])

# f = Figure()
# scatter([0, 1],[0, 1],[0, 1], plot_color=:red)
# Makie.plot(tq)
# Makie.plot!(tq)
# s = GLMakie.plot(tqp)

plotlyjs(reuse=true)
splot = Plots.surface(tqp, color=:blue, plot_caps=true, size=(1920,1080))
display(splot)
Plots.surface!(tqh, color=:red)
Plots.surface!(tqc, color=:green)
Plots.plot!([center center+3*ax][1,:], [center center+3*ax][2,:], [center center+3*ax][3,:], color=:blue, linewidth=4)
# current()

nray = 1000
rays = Array{Particle, 1}(undef, nray)

chits = zeros(1,nray)
phits = zeros(1,nray)
hhits = zeros(1,nray)

crays = zeros(3,nray)
prays = zeros(3,nray)
hrays = zeros(3,nray)

cnorm = zeros(3,nray)
pnorm = zeros(3,nray)
hnorm = zeros(3,nray)

for i = 1:nray
    dir = rand(3,) - 0.5*ones(3,)
    dir = dir/norm(dir)
    rays[i] = Particle(center + 1*ax, dir, 10)

    (cin, cout) = inout(rays[i], qc)
    (pin, pout) = inout(rays[i], qp)
    (hin, hout) = inout(rays[i], qh)
    chits[i] = max(cin, cout)
    phits[i] = max(pin, pout)
    hhits[i] = max(hin, hout)
    crays[:,i] = rays[i].r0 + chits[i]*rays[i].v
    prays[:,i] = rays[i].r0 + phits[i]*rays[i].v
    hrays[:,i] = rays[i].r0 + hhits[i]*rays[i].v
    cnorm[:,i] = -normal(qc, crays[:,i])
    pnorm[:,i] = normal(qp, prays[:,i])
    hnorm[:,i] = normal(qh, hrays[:,i])
end

# Plots.quiver(prays[1,:], prays[2,:], prays[3,:], quiver=(pnorm[1,:], pnorm[2,:], pnorm[3,:]))


Plots.scatter!([prays[1,:] prays[1,:]+pnorm[1,:]], [prays[2,:] prays[2,:]+pnorm[2,:]], [prays[3,:] prays[3,:]+pnorm[3,:]], color=:blue, markersize=1)

Plots.scatter!([hrays[1,:] hrays[1,:]+hnorm[1,:]], [hrays[2,:] hrays[2,:]+hnorm[2,:]], [hrays[3,:] hrays[3,:]+hnorm[3,:]], color=:red, markersize=1)

Plots.scatter!([crays[1,:] crays[1,:]+cnorm[1,:]], [crays[2,:] crays[2,:]+cnorm[2,:]], [crays[3,:] crays[3,:]+cnorm[3,:]], color=:green, markersize=1)


a = 9
Plots.xlims!((center[1]-a,center[1]+a))
Plots.ylims!((center[2]-a,center[2]+a))
Plots.zlims!((center[3]-a,center[3]+a))

display(splot)

