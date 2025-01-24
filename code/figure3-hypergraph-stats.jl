using Distributions, Clustering, NearestNeighbors, Random, Plots, LaTeXStrings, Colors
using Measures
using StatsBase
function hypergraph_edges(X,deg_list;radfunc=(dist,deg) -> dist/sqrt(deg))
    T = BallTree(X)
    edges = Vector{Int}[]
    for i in axes(X,2)
        # deg = min(ceil(Int,rand(degreedist)),n-1)
        deg = deg_list[i]
        idxs, dists = knn(T, X[:,i], deg+1)
        if deg > 1 
            maxdist = maximum(dists) 
            pts = @view X[:,idxs]
            rad = radfunc(maxdist,deg)
            clusters = dbscan(pts, rad).clusters
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
        # else
            # # only one vertex! 
            # push!(edges, [i,idxs[2]])
        end 
    end 
    return edges
end

function get_func(alpha) 
    if alpha <= 1
        return (d,deg) -> alpha*d/sqrt(deg)
    else
        return (d,deg) -> d/sqrt(deg) + (alpha-1)*(d-d/sqrt(deg))
    end
end 

Random.seed!(1234)

d = 2
n = 50000
alphas = range(0,2,length=25)
##### plot heatmap about hyperedge size distribution

X = rand(d,n)

deg_list = zeros(Int, n)
degreedist = LogNormal(log(3),1)
for i = 1:n
    deg_list[i] = min(ceil(Int,rand(degreedist)),n-1)
end

graphs = map(alpha -> hypergraph_edges(X, deg_list;radfunc=get_func(alpha)), alphas)

max_edge_size = mapreduce(e -> length(e), max, Iterators.flatten(graphs), init=2)

edge_cnt = zeros(Int, max_edge_size, 25)
for i = 1:25
    for e in graphs[i]
        edge_cnt[length(e), i] += 1
    end
end

edge_cnt_normalized = edge_cnt ./ sum(edge_cnt, dims=1)
### plot number of hyperedges vs alpha

function run_trial()
    X = rand(d,n)

    deg_list = zeros(Int, n)
    degreedist = LogNormal(log(3),1)
    for i = 1:n
        deg_list[i] = min(ceil(Int,rand(degreedist)),n-1)
    end

    graphs = map(alpha -> hypergraph_edges(X, deg_list;radfunc=get_func(alpha)), alphas)
    nedges = map(x->length(x), graphs)
    return nedges
end
trials = map(x->run_trial(), 1:25)  

nedges = reduce(hcat, trials) # now we have alpha as row index and trial as column index
# I want to take the quantiles over the columns... 
quantiles = [0.1,0.25,0.5,0.75,0.9]
linewidths = [0.5, 1, 2, 1, 0.5]
colors = [:lightgrey, :grey, :black, :grey, :lightgrey]

## PLOT 3 - alpha = 2 
Random.seed!(42)
X = rand(d,n)

deg_list = zeros(Int, n)
degreedist = LogNormal(log(3),1)
for i = 1:n
    deg_list[i] = min(ceil(Int,rand(degreedist)),n-1)
end

hedges = hypergraph_edges(X, deg_list;radfunc=get_func(2.0))


l = @layout [grid(1, 3, widths=[0.31, 0.38, 0.31])]

plt = Plots.plot(layout=l, margin=6*Plots.mm)
Plots.plot!(plt, size=(1200,300))

Plots.scatter!(plt,pairs(countmap(length.(hedges))),
                    yscale=:log10,
                    markerstrokewidth=0,
                    leg=false,
                    xscale=:log10,
                    xlabel="Number of Hyperedges",
                    ylabel="Hyperedge Size",
                    subplot=1,
                    # annotation=(90,5e3,L"\alpha=2")
                )
Plots.heatmap!(plt, log10.(edge_cnt_normalized), 
        xlabel=L"\alpha",
        ylabel="Hyperedge size",
        xticks=([1,7,13,19,25], [0.0, 0.5, 1.0, 1.5, 2.0]),
        # yticks=([1,2,3,4,5,6,7,8,9],[2,3,4,5,6,7,8,9,10]),
        yscale=:log10,
        # zscale=:log10, 
        color=:viridis,
        # color=:cividis,
        # color=:inferno,
        # color=:plasma,
        # color=:reds,
        subplot=2)


Plots.plot!(plt, xlabel=L"\alpha", ylabel="Number of hyperedges", subplot=3)
for (q, lw, c) in zip(quantiles, linewidths, colors)
    nedges_q = quantile.(eachrow(nedges), q)
    Plots.plot!(plt, alphas, nedges_q, label="", linewidth=lw, color=c, subplot=3)
end

#touch up margins
Plots.plot!(plt[2],
                left_margin=0Measures.mm,
                right_margin=0mm,
                bottom_margin=8mm
                )

display(plt)

plot!(plt,dpi=2000)
savefig(plt, "data/output/figures/final/hyperedge-degree-info.pdf")
savefig(plt, "data/output/figures/poster/hyperedge-degree-info.png")