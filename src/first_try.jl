using Contour
using Plots
using PolygonOps

include("example_data.jl")

dx = 0.1
xs = -10.0:dx:15.0
ys = -10.0:dx:20.0

phi = [two_gaussians_plus_negative_2d(x, y) for x in xs, y in ys]
p = plot()
# NB: it appears Contour expects the opposite shape to heatmap
heatmap!(p, xs, ys, phi')

for cl in levels(contours(xs, ys, phi, [-0.1, 0.0, 0.2]))
    lvl = level(cl)
    for line in lines(cl)
        xs_cnt, ys_cnt = coordinates(line)
        plot!(p, xs_cnt, ys_cnt)
    end
end
plot!(p)

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


function plot_area_levels(xs, ys, phi, levs)
    p_hm = heatmap(xs, ys, phi')
    areas = []
    area_levs = []
    for cl in levels(contours(xs, ys, phi, levs))
        lvl = level(cl)
        for line in lines(cl)
            pts = coordinates(line)
            plot!(p_hm, pts[1], pts[2], label=nothing, color=:black, linestyle=:dash)
            a = abs(area(coords2poly(pts)))
            push!(areas, a)
            push!(area_levs, lvl)
        end
    end
    plot(
        heatmap(xs, ys, phi'),
        p_hm,
        scatter(areas, levs, xlabel="area", ylabel="contour value", markersize=1, xlim=(0, 500^2)),
        layout=[1,1]
    )
end

plot_area_levels(xs, ys, phi, -0.9:0.02:0.0)
plot_area_levels(xs, ys, phi, -0.0:0.02:1.9)
plot_area_levels(xs, ys, phi, -0.9:0.02:1.9)

# test, area of a circle
dx = 0.1
xs = -2.0:dx:2.0
ys = -2.0:dx:2.0
phi = [sphere_2d(x, y, 1.0) for x in xs, y in ys]
pts = coordinates(lines(levels(contours(xs, ys, phi, [0.0]))[1])[1]);
isapprox(area(coords2poly(pts)), π; rtol=0.01)

## try with LES data

# raster = Raster("rico.no_shear_br0.05.qv.tn6.nc")
using NCDatasets

ds = NCDataset("rico.no_shear_br0.05.qv.tn6.nc")
zi = 6
zn = 36
phi = replace(ds[:qv][:,1,:,:][:,:], missing => 0.0)[zi:zn,:]
phi = replace(ds[:w][:,1,:,:][:,:], missing => 0.0)[zi:zn,:]
xs = replace(ds[:xt][:], missing => 0.0)
ys = replace(ds[:zt][:], missing => 0.0)[zi:zn]
heatmap(xs, ys, phi)

function find_area_levels(xs, ys, phi, levs)
    areas = []
    area_levs = []
    cnt_pts = []
    for cl in levels(contours(xs, ys, phi, levs))
        lvl = level(cl)
        for line in lines(cl)
            pts = coordinates(line)
            a = abs(area(coords2poly(pts)))
            push!(areas, a)
            push!(area_levs, lvl)
            push!(cnt_pts, pts)
        end
    end
    return (areas, area_levs, cnt_pts)
end

a_cnts, l_cnts, cnt_pts = find_area_levels(xs, ys, phi', 15.2:0.005:20.0);

plot(
    heatmap(xs ./ 1e3, ys, phi, ylabel="altitude [m]", xlabel="horizontal dist [km]"),
    scatter(sqrt.(a_cnts), l_cnts, xlabel="equiv rad. [m]", ylabel="qv contour [g/kg]", markersize=1, xlim=(0, 500)),
    layout=[1,2],
)

w = 2.0
p_pm = heatmap(xs ./ 1e3, ys, phi, ylabel="altitude [m]", xlabel="horizontal dist [km]"),
p_sc = scatter(sqrt.(a_cnts), l_cnts, xlabel="equiv rad. [m]", ylabel="qv contour [g/kg]", markersize=1, xlim=(0, 500))
for s in 200:50:550
    idx = s - w .< sqrt.(a_cnts) .< s + w;
    for (i, m) in enumerate(idx)
        if m
        plot!(p_hm, cnt_pts[i][1], cnt_pts[i][2], label=nothing, color=:black, linestyle=:dash)
        scatter!(p_sc, sqrt(a_cnts[i]), l_cnts[i], color=:red)
        end
    end
end
p = plot(size=(800, 100))
plot!(p)



function plot_area_levels(xs, ys, phi, levs)
    p_hm = heatmap(xs, ys, phi')
    areas = []
    area_levs = []
    for cl in levels(contours(xs, ys, phi, levs))
        lvl = level(cl)
        for line in lines(cl)
            pts = coordinates(line)
            plot!(p_hm, pts[1], pts[2], label=nothing, color=:black, linestyle=:dash)
            a = abs(area(coords2poly(pts)))
            push!(areas, a)
            push!(area_levs, lvl)
        end
    end
    plot(
        p_hm,
        scatter(sqrt.(areas), area_levs, xlabel="equiv rad. [m]", ylabel="qv contour [g/kg]", markersize=1, xlim=(0, 500)),
        layout=[1,1]
    )
end

plot_area_levels(xs, ys, phi', 15.0:0.01:19.0)

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