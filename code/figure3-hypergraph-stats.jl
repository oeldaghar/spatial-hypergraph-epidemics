using Distributed
# Distributed.addprocs(10)
@everywhere include("spatial-hypergraph.jl")
@everywhere using Distributions, Random
@everywhere using StatsBase
@everywhere using ProgressMeter
@everywhere using MatrixNetworks, Graphs, SparseArrays
@everywhere using SparseArrays

# Plotting 
using Plots, LaTeXStrings, Colors
using Measures

@everywhere function hedges_to_biadj(hedges,n)
    # hedges[hyperedge_id] = [nodes in hyperedge] 
    nedges = lastindex(hedges)
    
    nodes = Vector{Int}()
    edges = Vector{Int}()
    for (row,h) in enumerate(hedges)
        append!(nodes,h)
        append!(edges,row*ones(lastindex(h)))
    end
    return sparse(edges,nodes,ones(lastindex(nodes)),lastindex(hedges),n)
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

# @everywhere function hypergraph_edges(X,deg_list;radfunc=(dist,deg) -> dist/sqrt(deg))
#     T = BallTree(X)
#     edges = Vector{Int}[]
#     dim = size(X,1)
#     for i in axes(X,2)
#         # deg = min(ceil(Int,rand(degreedist)),n-1)
#         deg = deg_list[i]
#         idxs, dists = knn(T, X[:,i], deg+1)
#         if deg > 1 
#             maxdist = maximum(dists) 
#             pts = @view X[:,idxs]
#             rad = radfunc(maxdist,deg,dim)
#             clusters = dbscan(pts, rad).clusters
#             for c in clusters
#                 e = [i]
#                 for v in c.core_indices
#                     if idxs[v] != i
#                         push!(e, idxs[v])
#                     end
#                 end
#                 for v in c.boundary_indices
#                     if idxs[v] != i 
#                         push!(e, idxs[v])
#                     end
#                 end
#                 if length(e) > 1
#                     push!(edges, e)
#                 end
#             end
#         # else
#             # # only one vertex! 
#             # push!(edges, [i,idxs[2]])
#         end 
#     end 
#     return edges
# end

# @everywhere function get_func(alpha) 
#     if alpha <= 1
#         return (d,deg,dim) -> (alpha^(1/dim)*d)/(deg)^(1/(2*dim))
#     else
#         return (d,deg,dim) -> d/(deg)^(1/(2*dim)) + (alpha-1)*(d-d/(deg)^(1/(2*dim)))
#     end
# end 

@everywhere function get_edge_counts(n,d,alphas)
    # Random.seed!(1234)
    X = rand(d,n)
    
    deg_list = zeros(Int, n)
    degreedist = LogNormal(log(3),1)
    for i = 1:n
        deg_list[i] = min(ceil(Int,rand(degreedist)),n-1)
    end
    
    graphs = pmap(alpha -> hypergraph_edges(X, deg_list;radfunc=get_func(alpha)), alphas)
    
    max_edge_size = mapreduce(e -> length(e), max, Iterators.flatten(graphs), init=2)
    
    edge_cnt = zeros(Int, max_edge_size, lastindex(alphas))
    for i = 1:lastindex(alphas)
        for e in graphs[i]
            edge_cnt[length(e), i] += 1
        end
    end
    
    edge_cnt_normalized = edge_cnt ./ sum(edge_cnt, dims=1)
    return edge_cnt
end

@everywhere function get_edge_counts(graph)
    # Random.seed!(1234)    
    max_edge_size = maximum(length.(graph))
    
    edge_cnt = zeros(Int, max_edge_size)
    for e in graph
        edge_cnt[length(e)] += 1
    end
    
    edge_cnt_normalized = edge_cnt ./ sum(edge_cnt, dims=1)
    return edge_cnt
end

function sum_columns_of_matrix(matrix_of_vectors)
    max_length = maximum(length.(matrix_of_vectors))
    summed_matrix = zeros(size(matrix_of_vectors, 1), max_length)
    for (i, col) in enumerate(eachcol(matrix_of_vectors))
        for (j, vec) in enumerate(col)
            for k in 1:length(vec)
                summed_matrix[j, k] += vec[k]
            end
        end
    end
    return Matrix(summed_matrix')
end


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

@everywhere function run_trial(n,d,alphas,centers=1)
    X = rand(d,n)
    # X,_ = make_blobs(n,d,centers=centers,as_table = false)
    # X = Matrix(X')

    deg_list = zeros(Int, n)
    degreedist = LogNormal(log(3),1)
    for i = 1:n
        deg_list[i] = min(ceil(Int,rand(degreedist)),n-1)
    end

    graphs = pmap(alpha -> hypergraph_edges(X, deg_list;radfunc=get_func(alpha)), alphas)
    nedges = pmap(x->length(x), graphs)
    ntris = pmap(x->sum(1 for i in MatrixNetworks.triangles(pairwise_proj(x))), graphs)
    projected_degs = pmap(x->_project_graph_avgdegrees(x,scalefunc=y->1/y),graphs)
    edge_cnts = pmap(x->get_edge_counts(x),graphs)
    # global_cc_g = pmap(x->get_global_cc(x), graphs)
    return nedges,ntris,projected_degs,edge_cnts
end

function get_row_plotting_data(n=10000,d=2,alphas=range(0,2,25),ntrials=25,centers=1)
    # edge_cnts = get_edge_counts(n,d,range(0,2,25))

    trials = @showprogress map(x->run_trial(n,d,alphas,centers), 1:ntrials)  

    nedges = reduce(hcat, first.(trials)) # now we have alpha as row index and trial as column index
    ntris = reduce(hcat, map(x->x[2],trials))
    projected_degs = reduce(hcat, map(x->x[3],trials))
    # handle edge_cnts. comes in as a vector for each alpha and trial
    # matrix of alphas vs trial. vector in each component
    edge_cnts = reduce(hcat,map(x->x[4],trials))
    # aggreagte across trials #ncols 
    edge_cnts = sum_columns_of_matrix(edge_cnts)
    edge_cnts = edge_cnts./=sum(edge_cnts,dims=1)

    return [edge_cnts, nedges, ntris]
end 

## PLOTTING CODE 
function _custom_heatmap(pdata)
    # function to bin rows of the heatmap weight matrix
    function _bin_indices(bins::Vector, pdata::Vector)
        bin_dict = Dict{Float64, Vector{Int}}()
        for i in 1:length(bins)-1
            bin_dict[bins[i]] = findall(x -> bins[i] <= x < bins[i+1], pdata)
        end
        bin_dict[bins[end]] = findall(x -> x >= bins[end], pdata)
        return bin_dict
    end
    
    function _log_bin_ydata(pdata,ybins)
        binned_inds = _bin_indices(ybins,collect(1:lastindex(pdata,1)))
        new_mat = zeros(Float64,lastindex(ybins)-1,lastindex(pdata,2))
        for (newrow,key) in enumerate(ybins[1:end-1])
            old_rows = binned_inds[key]
            new_mat[newrow,:] = sum(pdata[old_rows,:],dims=1)
        end
        return new_mat
    end
    
    max_hsize = lastindex(pdata,1)
    ybins = (10.0).^(range(1,log10(max_hsize+10),15))
    ybins = vcat(1:9,ybins)
    ybins = sort(unique(ybins))

    new_mat = _log_bin_ydata(pdata,ybins)

    # xrange, yrange, data_matrix, yaxis scale
    f = Plots.heatmap(1:lastindex(new_mat,2),ybins[1:end-1],log10.(new_mat),
                yscale=:log10,
                color=:viridis,
                clims=(-5,0),
                )
    return f#,new_mat,ybins 
end

function make_fig(data)
    # alphas = range(0,2,15)
    # n = 5000
    # println("Simulating data...")
    # row1 = get_row_plotting_data(n,2,alphas)
    # row2 = get_row_plotting_data(n,5,alphas)
    # row3 = get_row_plotting_data(n,10,alphas)
    # data = [row1,row2,row3]

    # handle heatmaps
    # make same size 
    max_hyperedge_size = maximum(map(x->size(x[1],1),data))
    for row in data
        heatmap_data = row[1]
        rows_to_pad = max_hyperedge_size-size(heatmap_data,1)
        row[1] = vcat(heatmap_data,zeros(rows_to_pad,size(heatmap_data,2)))
    end

    # make individual figures 
    # heatmaps 
    col1_figs = []
    for pdata in first.(data)
        f = _custom_heatmap(pdata)
        Plots.plot!(f,
                xlabel=L"\alpha",
                ylabel="Hyperedge size",
                xticks=([1,7,13,19,25], [0.0, 0.5, 1.0, 1.5, 2.0]),
                yscale=:log10,
                color=:viridis,
                # clims=(-6.0,0),
                framestyle=:box,
                thickness_scaling=1.2,
                guidefontsize=14,
                tickfontsize=12)
        # f = Plots.heatmap(log10.(pdata), 
        #         xlabel=L"\alpha",
        #         ylabel="Hyperedge size",
        #         xticks=([1,7,13,19,25], [0.0, 0.5, 1.0, 1.5, 2.0]),
        #         yscale=:log10,
        #         color=:viridis,
        #         clims=(-5.0,0),
        #         framestyle=:box,
        #         thickness_scaling=1.2,
        #         guidefontsize=14,
        #         tickfontsize=12)
        push!(col1_figs,f)
    end

    # col2 and col3 
    quantiles = [0.1,0.25,0.5,0.75,0.9]
    linewidths = [0.5, 1, 2, 1, 0.5]
    colors = [:lightgrey, :grey, :black, :grey, :lightgrey]

    # total edges 
    col2_figs = []
    for pdata in map(x->x[2],data)
        f = Plots.plot(xlabel=L"\alpha", ylabel="Total hyperedges", 
                    framestyle=:box,
                    thickness_scaling=1.2,
                    guidefontsize=14,
                    tickfontsize=12,
                    tickdirection=:out)
        for (q, lw, c) in zip(quantiles, linewidths, colors)
            nedges_q = quantile.(eachrow(pdata), q)
            Plots.plot!(f, alphas, nedges_q, label="", linewidth=lw, color=c,    
                yscale=:log10)
        end
        push!(col2_figs,f)
    end
    # align ylims 
    col_ylims = Plots.ylims(Plots.plot(col2_figs...,layout=(3,1),link=:all))
    for f in col2_figs
        Plots.plot!(f,ylims=col_ylims)
    end

    col3_figs = []
    for pdata in map(x->x[3],data)
        f = Plots.plot(xlabel=L"\alpha", ylabel="Total Triangles", 
                    framestyle=:box,
                    thickness_scaling=1.2,
                    guidefontsize=14,
                    tickfontsize=12,
                    tickdirection=:out)
        for (q, lw, c) in zip(quantiles, linewidths, colors)
            ntris_q = quantile.(eachrow(pdata), q)
            Plots.plot!(f, alphas, ntris_q, label="", linewidth=lw, color=c,    
                yscale=:log10)
        end
        push!(col3_figs,f)
    end
    # align ylims 
    col_ylims = Plots.ylims(Plots.plot(col3_figs...,layout=(3,1),link=:all))
    for f in col3_figs
        Plots.plot!(f,ylims=col_ylims)
    end


    # put them all together 
    figs =[]
    for tup in zip(col1_figs,col2_figs,col3_figs)
        push!(figs,tup[1])
        push!(figs,tup[2])
        push!(figs,tup[3])
    end

    l = @layout [Plots.grid(3, 3, widths=[0.42, 0.29, 0.29])]
    plt = Plots.plot(figs...,layout=l, 
                margin=0*Plots.mm, size=(1200,1100))
    Plots.plot!(plt,top_margin = 10mm,bottom_margin=-2mm)
    Plots.plot!(plt[2],title="n=$n  d=2",titlefontsize = 24,
                    top_margin=-3Measures.mm)
    Plots.plot!(plt[5],title="n=$n  d=5",titlefontsize = 24,
                    top_margin=-3Measures.mm)
    Plots.plot!(plt[8],title="n=$n  d=10",titlefontsize = 24,
                    top_margin=-3Measures.mm)
    return plt 
end

alphas = range(0,2,25)
n = 10000
centers = 1
println("Simulating data...")
row1 = get_row_plotting_data(n,2,alphas,centers)
row2 = get_row_plotting_data(n,5,alphas,centers)
row3 = get_row_plotting_data(n,10,alphas,centers)
data = [row1,row2,row3]

plt = make_fig(data)

Plots.savefig(plt,"data/output/figures/final/hypergraph-stats.pdf")
# Plots.plot!(plt,dpi=1000)
# Plots.savefig(plt,"data/output/figures/final/hypergraph-stats.png")
