using Distributions
using Statistics
using NearestNeighbors
using Clustering
using LazySets
using Plots
using Random

using Measures

function circle(x, y, r)
    theta = LinRange(0, 2*pi, 500)
    x .+ r * sin.(theta), y .+ r * cos.(theta)
end

plt = Plots.plot(xlims=(0,1),
           ylims=(0,1),
           aspect_ratio=:equal,
           framestyle=:none,
           background_color_inside=:gray95,
           legend=false,
           ticks=false,
           layout=(1,3))

Plots.plot!(plt, size=(800,200))

xy = [[0.5,0.5] [0.38,0.45] [0.44,0.36] [0.71,0.71] [0.64,0.78] [0.46,0.72] [0.56,0.67] [0.82,0.25] [0.72,0.15]]

annotation_offsets = [(-0.05,0.05),(-0.06,-0.01),(-0.05,-0.05),(0.05,0.05),(0.01,0.07),(-0.05,0.05),(-0.01,0.07),(0.01,0.07),(-0.06,-0.02)]
annotation_labels = [Plots.text(i, 11) for i=1:lastindex(annotation_offsets)]

n = size(xy, 2)

T = BallTree(xy)

i = 1
deg = n-1
idxs, dists = knn(T, xy[:,i], deg+1)
pts = @view xy[:,idxs]

rad = 0.15

function make_annotations(f,krange=range(1,n))
    for k in krange 
        xoffset,yoffset = annotation_offsets[k][1],annotation_offsets[k][2]
        annotate!(f,xy[1,k]+xoffset,xy[2,k]+yoffset,annotation_labels[k])
    end
end

function plot_points(f)
    Plots.scatter!(f, 
        xy[1,2:end], 
        xy[2,2:end], 
        color=1,
        markerstrokewidth=0,
        markersize=5)
    Plots.scatter!(f, 
        xy[1,1:1], 
        xy[2,1:1], 
        color=:red,
        markerstrokewidth=0,
        markersize=5)
end

#figure 1 points with circle but no edges 
plt1 = Plots.plot(size = (200,200),
            xlims=(0,1),
           ylims=(0,1),
           aspect_ratio=:equal,
           framestyle=:none,
           background_color_inside=:gray95,
           legend=false,
           ticks=false
           )
make_annotations(plt1)
Plots.plot!(plt1, 
    circle(xy[1,i], xy[2,i], 0.45), 
    c=:red)
plot_points(plt1)
plt1

###FIGURE 2 - cluster neighbors
plt2 = deepcopy(plt1)
#cluster neighbors 
clusters = dbscan(pts, rad)
edges = Vector{Int}[]
for c in clusters
    e = [i]
    for v in c.core_indices
        if idxs[v] != i
            push!(e, idxs[v])
        end
    end
    for v in c.boundary_indices
        if idxs[v] != i 
            push!(e, idxs[v])
        end
    end
    if length(e) > 1
        push!(edges, e)
    end
end
#highlight edges 
for edge in edges
    e = [x for x in edge if x!=1]
    if length(e) >= 2
        local b = [Ball2(xy[:, v], 0.02) for v in e]
        local c = ConvexHullArray(b)
        Plots.plot!(plt2, c, 1e-3, alpha=0.1, c=:blue)
    else
        Plots.plot!(plt2, xy[1,e], xy[2,e], c=:blue)
    end
end
plt2

### FIGURE 3
plt3 = Plots.plot(size = (200,200),
            xlims=(0,1),
           ylims=(0,1),
           aspect_ratio=:equal,
           framestyle=:none,
           background_color_inside=:gray95,
           legend=false,
           ticks=false
           )
Plots.plot!(plt3, 
    circle(xy[1,i], xy[2,i], 0.45), 
    c=:red)

for edge in edges
    e = edge
    if length(e) >= 2
        local b = [Ball2(xy[:, v], 0.02) for v in e]
        local c = ConvexHullArray(b)
        Plots.plot!(plt3, c, 1e-3, alpha=0.1, c=:blue)
    else
        Plots.plot!(plt3, xy[1,e], xy[2,e], c=:blue)
    end
end

make_annotations(plt3)
plot_points(plt3)

## putting it together
plt = plot(plt1,plt2,plt3,layout=(1,3),size=(600,200))
plot!(plt,margins = -1Measures.mm)

savefig(plt, "data/output/figures/final/fig0.pdf")
# savefig(plt, "data/output/figures/poster/fig0.svg")