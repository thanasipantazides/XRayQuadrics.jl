using XRayTrace

n = [0;0;1]
c = [1;1;1]

p = Plane(c,n)
q = Quadric(p)


ray = XRayTrace.Particle([0;0;10],[0;-1;-1],10e3)

ts = XRayTrace.inout(ray, q)