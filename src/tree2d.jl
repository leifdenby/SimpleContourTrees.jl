# Relevant discussion:
# - https://discourse.julialang.org/t/most-effective-way-to-check-if-neighboring-values-in-a-matrix-are-the-same/86291/12
using Plots

dx = 0.1
xs = 0.0:dx:10.0
ys = 0.0:dx:10.0

gauss_2d(x, y, σ; x0=0.0, y0=0.0) = exp(-((x-x0)^2 + (y-y0)^2)/σ^2)
sphere_2d(x, y, r0; x0=0.0, y0=0.0) = (x-x0)^2 + (y-y0)^2 - r0^2

function two_gaussians_2d(x, y)
    w = 3.0
    gauss_2d(x, y, w, x0=2.0, y0=2.0) + gauss_2d(x, y, w, x0=8.0, y0=7.0)
end

function two_gaussians_plus_negative_2d(x, y)
    w = 3.0
    two_gaussians_2d(x, y) - gauss_2d(x, y, w, x0=2.0, y0=7.0)
end

phi = [two_gaussians_2d(x, y) for x in xs, y in ys]

heatmap(xs, ys, phi)


NeighborIdxs2D = [
    CartesianIndex(i, j)
    for (i, j) in [(-1, 0), (1, 0), (0, -1), (0, 1)]
]

NeighborIdxs3D = [
    CartesianIndex(i, j, k)
    for (i, j, k) in [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)]
]

function neighbor_offsets(N)
    if N == 2
        return NeighborIdxs2D
    elseif N == 3
        return NeighborIdxs3D
    else
        throw("$N not implemented")
    end
end

"""
Find extrema by looking at all neighbors in 2D or 3D array.

Returns array of same shape is input with [-1, 0, 1] indicating minima, not
extrema and maxima respectively
"""
function find_extrema(ϕ::AbstractArray{T,N}) where {T,N}
    out = zero(ϕ)
    I_inset = CartesianIndices(tuple([2:size(ϕ, i)-2 for i in 1:ndims(ϕ)]...))
    Ineighbors = neighbor_offsets(N)
    for I in I_inset
        extrema_kind = zero(T)
        ϕ_local = ϕ[I]
        for J in Ineighbors
            sign_local = sign(ϕ_local - ϕ[I+J]) 
            if extrema_kind == 0.0
                extrema_kind = sign_local
            elseif extrema_kind != sign_local
                extrema_kind = zero(T)
                break
            end
        end
        out[I] = extrema_kind
    end
    out 
end

heatmap(find_extrema(phi))


function assign_to_extrema(phi::AbstractArray{T,N}, extrema) where {T, N}
    out = zero(extrema)
    extrema_list = []
    for i in UnitRange(eachindex(phi))
        if extrema[i] != 0
            out[i] = i
            push!(extrema_list, CartesianIndices(phi)[i])
        end
    end
    
    cart2idx = reshape(eachindex(phi), size(phi))

    Ineighbors = neighbor_offsets(N)
    for I_ext in extrema_list
        kind = extrema[I_ext]
        new_points = [I_ext]
            
        n = 0
        while length(new_points) > 0
            next_points = copy(new_points)
            new_points = []
            for I_c in next_points
                n += 1
                for J in Ineighbors
                    I = I_c + J
                    try
                        if out[I] != 0
                            continue
                        elseif sign(phi[I_c] - phi[I]) == kind
                            out[I] = cart2idx[I_ext]
                            push!(new_points, I)
                        end
                    catch
                    end
                end
            end
            if n > 6000
                break
            end
        end
    end
    out
end

extrema = find_extrema(phi)
extrema = replace(extrema, -1 => 0)
plot(
    heatmap(phi),
    heatmap(extrema),
    heatmap(assign_to_extrema(phi, extrema))
)
