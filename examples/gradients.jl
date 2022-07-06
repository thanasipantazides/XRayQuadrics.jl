using XRayQuadrics
using LinearAlgebra, StaticArrays
using Profile

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

randv = rand(3)
randv = randv/norm(randv)

r = SVector{3}(R*(randv×axis)/norm(randv×axis))

@time nd = normal(qc, r)
@time ng = normal_grad(qc, r)

println(nd'*ng ≈ 1.0)