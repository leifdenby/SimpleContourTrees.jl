using Contour
using Plots
using PolygonOps

dx = 0.01
xs = 0.0:dx:4.0
ys = 0.0:dx:4.0

gauss_2d(x, y, σ; x0=0.0, y0=0.0) = exp(-((x-x0)^2 + (y-y0)^2)/σ)
sphere_2d(x, y, r0; x0=0.0, y0=0.0) = (x-x0)^2 + (y-y0)^2 - r0^2

phi = [gauss_2d(x, y, 0.5, x0=0.5, y0=0.5) for x in xs, y in ys]
phi = [sphere_2d(x, y, 1.0, x0=2.0, y0=2.0) for x in xs, y in ys]

p = plot()
heatmap!(p, xs, ys, phi)

for cl in levels(contours(xs, ys, phi, [0.0,]))
    lvl = level(cl)
    for line in lines(cl)
        xs_cnt, ys_cnt = coordinates(line)
        plot!(p, xs_cnt, ys_cnt)
    end
end
plot!(p)

pts = coordinates(lines(levels(contours(xs, ys, phi, [0.0]))[1])[1]);

coords2poly(pts) = collect(zip(pts[1], pts[2]))

function area(poly::T) where T
    a = zero(eltype(eltype(T)))
    @inbounds for i in 1:length(poly)-1
        p1 = poly[i]
        p2 = poly[i+1]
        a += p1[1]*p2[2]-p2[1]*p1[2]
    end
    return a/2
end

isapprox(area(coords2poly(pts)), π; rtol=0.01)

## scratch space below

area(collect(zip(pts[1], pts[2])))

pts = [
    0 0 1 1
    0 1 1 0
]
pts = (pts[1,:], pts[2,:])

pts'

plot(pts[1,:], pts[2,:])

area(collect(zip(pts[1], pts[2])))

pts[1,:]


pts[1,:]


area(stack(pts))