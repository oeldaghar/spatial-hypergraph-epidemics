using Distributions
using Statistics
using NearestNeighbors
using Clustering
using LazySets
using Plots
using Random

function circle(x, y, r)
    theta = LinRange(0, 2*pi, 500)
    x .+ r * sin.(theta), y .+ r * cos.(theta)
end

plt = plot(xlims=(0,1),
           ylims=(0,1),
           aspect_ratio=:equal,
           framestyle=:none,
           background_color_inside=:gray95,
           legend=false,
           ticks=false,
           layout=(1,4))

plot!(plt, size=(800,200))

xy = [[0.5,0.5] [0.38,0.45] [0.44,0.36] [0.71,0.71] [0.64,0.78] [0.46,0.72] [0.56,0.67] [0.82,0.25] [0.72,0.15]]

n = size(xy, 2)

T = BallTree(xy)

i = 1
deg = n-1
idxs, dists = knn(T, xy[:,i], deg+1)
pts = @view xy[:,idxs]

function plot_subfig(subfig, rad)

    scatter!(plt, 
             xy[1,:], 
             xy[2,:], 
            #  series_annotations=text.(1:n, :bottom, 9),
             subplot=subfig)

    annotate!(plt, xy[1,1]-0.05, xy[2,1]+0.05, text(1, 11))
    annotate!(plt, xy[1,2]-0.06, xy[2,2]-0.01, text(2, 11))
    annotate!(plt, xy[1,3]-0.05, xy[2,3]-0.05, text(3, 11))
    annotate!(plt, xy[1,4]+0.05, xy[2,4]+0.05, text(4, 11))
    annotate!(plt, xy[1,5]+0.01, xy[2,5]+0.07, text(5, 11))
    annotate!(plt, xy[1,6]-0.05, xy[2,6]+0.05, text(6, 11))
    annotate!(plt, xy[1,7]-0.01, xy[2,7]+0.07, text(7, 11))
    annotate!(plt, xy[1,8]+0.01, xy[2,8]+0.07, text(8, 11))
    annotate!(plt, xy[1,9]-0.06, xy[2,9]-0.02, text(9, 11))

    clusters = dbscan(pts, rad).clusters
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

    plot!(plt, 
          circle(xy[1,i], xy[2,i], 0.45), 
          c=:red,
          subplot=subfig)

    for e in edges
        if length(e) > 2
            local b = [Ball2(xy[:, v], 0.02) for v in e]
            local c = ConvexHullArray(b)
            plot!(plt, c, 1e-3, alpha=0.1, subplot=subfig, c=:blue)
        else
            plot!(plt, xy[1,e], xy[2,e], subplot=subfig, c=:blue)
        end
    end
end

plot_subfig(1, 0.01)
plot_subfig(2, 0.12)
plot_subfig(3, 0.15)
plot_subfig(4, 0.4)

display(plt)

savefig(plt, "fig1.pdf")






