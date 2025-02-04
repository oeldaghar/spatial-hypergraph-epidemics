using Distributed
# Distributed.addprocs(10)

@everywhere using Distributions, Clustering, NearestNeighbors, Random, Plots, LaTeXStrings, Colors
@everywhere using StatsBase
@everywhere using ProgressMeter
@everywhere using MatrixNetworks, Graphs, SparseArrays
using Measures

@everywhere using SparseArrays
@everywhere function hedges_to_biadj(hedges,n)
    # hedges[hyperedge_id] = [nodes in hyperedge] 
    nedges = lastindex(hedges)
    
    nodes = Vector{Int}()
    edges = Vector{Int}()
    for (col,h) in enumerate(hedges)
        append!(nodes,h)
        append!(edges,col*ones(lastindex(h)))
    end
    return sparse(nodes,edges,ones(lastindex(nodes)),n,lastindex(hedges))
end

@everywhere function pairwise_proj(hedges)
    ei = Vector{Int}()
    ej = Vector{Int}()
    for h in hedges
        for node1 in h 
            for node2 in h 
                push!(ei,node1)
                push!(ej,node2)
            end
        end
    end
    nnodes = maximum([maximum(ei),maximum(ej)])
    A = sparse(ei,ej,ones(lastindex(ei)),nnodes,nnodes)
    for i=1:nnodes
        A[i,i]=0
    end
    dropzeros!(A)
    return A
end

@everywhere function hypergraph_edges(X,deg_list;radfunc=(dist,deg) -> dist/sqrt(deg))
    T = BallTree(X)
    edges = Vector{Int}[]
    dim = size(X,1)
    for i in axes(X,2)
        # deg = min(ceil(Int,rand(degreedist)),n-1)
        deg = deg_list[i]
        idxs, dists = knn(T, X[:,i], deg+1)
        if deg > 1 
            maxdist = maximum(dists) 
            pts = @view X[:,idxs]
            rad = radfunc(maxdist,deg,dim)
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

@everywhere function get_func(alpha) 
    if alpha <= 1
        return (d,deg,dim) -> (alpha^(1/dim)*d)/(deg)^(1/(2*dim))
        # return (d,deg,dim) -> alpha*d/(deg)^(1/(2))
    else
        return (d,deg,dim) -> d/(deg)^(1/(2*dim)) + (alpha-1)*(d-d/(deg)^(1/(2*dim)))
        # return (d,deg,dim) -> d/(deg)^(1/(2)) + (alpha-1)*(d-d/(deg)^(1/(2)))
    end
end 

Random.seed!(1234)

@everywhere d = 2
@everywhere n = 50000
@everywhere alphas = range(0,2,length=25)
##### plot heatmap about hyperedge size distribution

X = rand(d,n)
@everywhere X = $X

deg_list = zeros(Int, n)
degreedist = LogNormal(log(3),1)
for i = 1:n
    deg_list[i] = min(ceil(Int,rand(degreedist)),n-1)
end

graphs = pmap(alpha -> hypergraph_edges(X, deg_list;radfunc=get_func(alpha)), alphas)

max_edge_size = mapreduce(e -> length(e), max, Iterators.flatten(graphs), init=2)

edge_cnt = zeros(Int, max_edge_size, 25)
for i = 1:25
    for e in graphs[i]
        edge_cnt[length(e), i] += 1
    end
end

edge_cnt_normalized = edge_cnt ./ sum(edge_cnt, dims=1)
### plot number of hyperedges vs alpha

@everywhere function get_global_cc(hedges)
    g = SimpleGraph(pairwise_proj(hedges))
    return global_clustering_coefficient(g)
end 
@everywhere function _project_graph_avgdegrees(edges;scalefunc=x->1/x)
    sumdegs = 0.0
    for e in edges
        sumdegs += length(e)*(length(e)-1)*scalefunc(length(e))
    end
    return sumdegs
end

@everywhere function run_trial()
    X = rand(d,n)

    deg_list = zeros(Int, n)
    degreedist = LogNormal(log(3),1)
    for i = 1:n
        deg_list[i] = min(ceil(Int,rand(degreedist)),n-1)
    end

    graphs = pmap(alpha -> hypergraph_edges(X, deg_list;radfunc=get_func(alpha)), alphas)
    nedges = pmap(x->length(x), graphs)
    ntris = pmap(x->sum(1 for i in MatrixNetworks.triangles(pairwise_proj(x))), graphs)
    projected_degs = pmap(x->_project_graph_avgdegrees(x,scalefunc=y->1/y),graphs)
    # global_cc_g = pmap(x->get_global_cc(x), graphs)
    return nedges,ntris,projected_degs
end

trials = @showprogress map(x->run_trial(), 1:10)  

nedges = reduce(hcat, first.(trials)) # now we have alpha as row index and trial as column index
ntris = reduce(hcat, map(x->x[2],trials))
projected_degs = reduce(hcat, map(x->x[3],trials))

# I want to take the quantiles over the columns... 
quantiles = [0.1,0.25,0.5,0.75,0.9]
linewidths = [0.5, 1, 2, 1, 0.5]
colors = [:lightgrey, :grey, :black, :grey, :lightgrey]

## PLOT 3 - alpha = 2 
l = @layout [Plots.grid(1, 3, widths=[0.42, 0.29, 0.29])]

plt = Plots.plot(layout=l, margin=6*Plots.mm)
Plots.plot!(plt, size=(1200,300))

Plots.heatmap!(plt, log10.(edge_cnt_normalized), 
        xlabel=L"\alpha",
        ylabel="Hyperedge size",
        xticks=([1,7,13,19,25], [0.0, 0.5, 1.0, 1.5, 2.0]),
        yscale=:log10,
        color=:viridis,
        framestyle=:box,
        thickness_scaling=1.1,
        guidefontsize=14,
        tickfontsize=12,
        subplot=1)

Plots.plot!(plt, xlabel=L"\alpha", ylabel="Total hyperedges", 
            framestyle=:box,
            thickness_scaling=1.1,
            guidefontsize=14,
            tickfontsize=12,
            tickdirection=:out,
            subplot=2)
for (q, lw, c) in zip(quantiles, linewidths, colors)
    nedges_q = quantile.(eachrow(nedges), q)
    Plots.plot!(plt, alphas, nedges_q, label="", linewidth=lw, color=c, subplot=2,    
        yscale=:log10)
end

Plots.plot!(plt, xlabel=L"\alpha", ylabel="Projected Volume", 
            framestyle=:box,
            thickness_scaling=1.1,
            guidefontsize=14,
            tickfontsize=12,
            tickdirection=:out,
            subplot=3)
for (q, lw, c) in zip(quantiles, linewidths, colors)
    projected_degs_q = quantile.(eachrow(projected_degs), q)
    Plots.plot!(plt, alphas, projected_degs_q, label="", linewidth=lw, color=c, subplot=3)
end

#touch up margins
Plots.plot!(plt[2],
                left_margin=-3Measures.mm,
                right_margin=0mm,
                bottom_margin=8mm,
                )
Plots.plot!(plt[3],right_margin=1mm,)

Plots.plot!(plt,top_margin = 2mm)
display(plt)

plot!(plt,dpi=2000)
savefig(plt, "data/output/figures/final/hyperedge-stats-$n-$d-newalpha.pdf")
savefig(plt, "data/output/figures/final/hyperedge-stats-$n-$d-newalpha.png")