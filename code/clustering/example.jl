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

plt = Plots.plot(xlims=(0,1),
           ylims=(0,1),
           aspect_ratio=:equal,
           framestyle=:none,
           background_color_inside=:gray95,
           legend=false,
           ticks=false,
           layout=(3,2))

Plots.plot!(plt, size=(400,600))

xy = [[0.33,0.5] [0.6,0.5] [0.15,0.55] [0.12,0.4] [0.2,0.25] [0.83,0.7] [0.87,0.6] [0.76,0.34]]

v_annotation = (xy[1,1], xy[2,1]+0.07, Plots.text('v', 11))
u_annotation = (xy[1,2]-0.04, xy[2,2]-0.04, Plots.text('u', 11))


n = size(xy, 2)

T = BallTree(xy)

# the 3th input parameter i is the node id
# this function plot the hyperedge formed around node i (the i-th one in the given xy) whose degree is given by deg
# if alpha < 1, modify the function of rad
function plot_subfig(subfig, deg, i, alpha)

    idxs, dists = knn(T, xy[:,i], deg+1)
    maxdist = maximum(dists)
    pts = @view xy[:,idxs]
    rad = maxdist/sqrt(deg) + (alpha-1)*(maxdist-maxdist/sqrt(deg)) # alpha > 1
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

    Plots.plot!(plt, 
        circle(xy[1,i], xy[2,i], maxdist+0.015), 
        c=:red,
        subplot=subfig)

    for e in edges
        if length(e) > 2
            local b = [Ball2(xy[:, v], 0.02) for v in e]
            local c = ConvexHullArray(b)
            Plots.plot!(plt, c, 1e-3, alpha=0.1, subplot=subfig, c=:blue)
        else
            Plots.plot!(plt, xy[1,e], xy[2,e], subplot=subfig, c=:blue)
        end
    end

    Plots.scatter!(plt, 
            xy[1,:], 
            xy[2,:], 
            subplot=subfig,
            color=1,
            markerstrokewidth=0)

    annotate!(plt, v_annotation..., subplot=subfig)
    annotate!(plt, u_annotation..., subplot=subfig)
end

plot_subfig(1, 4, 1, 1.9)
plot_subfig(1, 4, 2, 1.9)
plot_subfig(2, 4, 1, 2)
plot_subfig(2, 4, 2, 2)

# plot subfigure 3
subfig = 3

# plot graph edges
Plots.plot!(plt, xy[1,[1,2]], xy[2,[1,2]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[1,3]], xy[2,[1,3]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[1,4]], xy[2,[1,4]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[1,5]], xy[2,[1,5]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[3,4]], xy[2,[3,4]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[3,5]], xy[2,[3,5]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[4,5]], xy[2,[4,5]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[1,6]], xy[2,[1,6]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[1,7]], xy[2,[1,7]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[1,8]], xy[2,[1,8]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[2,6]], xy[2,[2,6]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[2,7]], xy[2,[2,7]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[2,8]], xy[2,[2,8]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[6,7]], xy[2,[6,7]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[6,8]], xy[2,[6,8]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[7,8]], xy[2,[7,8]], subplot=subfig, c=:blue)

# plot graph nodes
Plots.scatter!(plt, 
        xy[1,:], 
        xy[2,:], 
        subplot=subfig,
        color=1,
        markerstrokewidth=0)

annotate!(plt, v_annotation..., subplot=subfig)
annotate!(plt, u_annotation..., subplot=subfig)

# plot subfigure 4
subfig = 4

Plots.plot!(plt, xy[1,[1,2]], xy[2,[1,2]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[1,3]], xy[2,[1,3]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[1,4]], xy[2,[1,4]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[1,5]], xy[2,[1,5]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[3,4]], xy[2,[3,4]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[3,5]], xy[2,[3,5]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[4,5]], xy[2,[4,5]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[1,6]], xy[2,[1,6]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[1,7]], xy[2,[1,7]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[1,8]], xy[2,[1,8]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[2,6]], xy[2,[2,6]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[2,7]], xy[2,[2,7]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[2,8]], xy[2,[2,8]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[6,7]], xy[2,[6,7]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[6,8]], xy[2,[6,8]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[7,8]], xy[2,[7,8]], subplot=subfig, c=:blue)
Plots.plot!(plt, xy[1,[2,3]], xy[2,[2,3]], subplot=subfig, c=:black)
Plots.plot!(plt, xy[1,[2,4]], xy[2,[2,4]], subplot=subfig, c=:black)
Plots.plot!(plt, xy[1,[2,5]], xy[2,[2,5]], subplot=subfig, c=:black)

Plots.scatter!(plt, 
        xy[1,:], 
        xy[2,:], 
        subplot=subfig,
        color=1,
        markerstrokewidth=0)

annotate!(plt, v_annotation..., subplot=subfig)
annotate!(plt, u_annotation..., subplot=subfig)



# now highlight the wedges of u 
subfig = 5

# plot graph edges
Plots.plot!(plt, xy[1,[1,2]], xy[2,[1,2]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[1,3]], xy[2,[1,3]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[1,4]], xy[2,[1,4]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[1,5]], xy[2,[1,5]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[3,4]], xy[2,[3,4]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[3,5]], xy[2,[3,5]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[4,5]], xy[2,[4,5]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[1,6]], xy[2,[1,6]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[1,7]], xy[2,[1,7]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[1,8]], xy[2,[1,8]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[2,6]], xy[2,[2,6]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[2,7]], xy[2,[2,7]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[2,8]], xy[2,[2,8]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[6,7]], xy[2,[6,7]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[6,8]], xy[2,[6,8]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[7,8]], xy[2,[7,8]], subplot=subfig, c=:grey)

# plot graph nodes
Plots.scatter!(plt, 
        xy[1,:], 
        xy[2,:], 
        subplot=subfig,
        color=1,
        markerstrokewidth=0)

Plots.scatter!(plt, 
    xy[1,[2]], 
    xy[2,[2]], 
    subplot=subfig,
    color=:black,
    markersize=6,
    markerstrokewidth=0)

annotate!(plt, v_annotation..., subplot=subfig)
annotate!(plt, u_annotation..., subplot=subfig)


subfig = 6

Plots.plot!(plt, xy[1,[1,2]], xy[2,[1,2]], subplot=subfig, c=:grey)

Plots.plot!(plt, xy[1,[1,3]], xy[2,[1,3]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[1,4]], xy[2,[1,4]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[1,5]], xy[2,[1,5]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[3,4]], xy[2,[3,4]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[3,5]], xy[2,[3,5]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[4,5]], xy[2,[4,5]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[1,6]], xy[2,[1,6]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[1,7]], xy[2,[1,7]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[1,8]], xy[2,[1,8]], subplot=subfig, c=:grey)

Plots.plot!(plt, xy[1,[2,6]], xy[2,[2,6]], subplot=subfig, c=2, linewidth=3)
Plots.plot!(plt, xy[1,[2,7]], xy[2,[2,7]], subplot=subfig, c=2, linewidth=3)
Plots.plot!(plt, xy[1,[2,8]], xy[2,[2,8]], subplot=subfig, c=2, linewidth=3)

Plots.plot!(plt, xy[1,[6,7]], xy[2,[6,7]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[6,8]], xy[2,[6,8]], subplot=subfig, c=:grey)
Plots.plot!(plt, xy[1,[7,8]], xy[2,[7,8]], subplot=subfig, c=:grey)

Plots.plot!(plt, xy[1,[2,3]], xy[2,[2,3]], subplot=subfig, c=3, linewidth=3)
Plots.plot!(plt, xy[1,[2,4]], xy[2,[2,4]], subplot=subfig, c=3, linewidth=3)
Plots.plot!(plt, xy[1,[2,5]], xy[2,[2,5]], subplot=subfig, c=3, linewidth=3)

Plots.scatter!(plt, 
        xy[1,[3,4,5]], 
        xy[2,[3,4,5]], 
        subplot=subfig,
        color=3,
        markerstrokewidth=0)

Plots.scatter!(plt, 
    xy[1,[1]], 
    xy[2,[1]], 
    subplot=subfig,
    color=1,
    markerstrokewidth=0)
    
Plots.scatter!(plt, 
    xy[1,[2]], 
    xy[2,[2]], 
    subplot=subfig,
    color=:black,
    markersize=6,
    markerstrokewidth=0)


Plots.scatter!(plt, 
    xy[1,[6,7,8]], 
    xy[2,[6,7,8]], 
    subplot=subfig,
    color=2,
    markerstrokewidth=0)


annotate!(plt, v_annotation..., subplot=subfig)
annotate!(plt, u_annotation..., subplot=subfig)


display(plt)
# TODO fix up the margins