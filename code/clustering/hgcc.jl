using Distributed 
# addprocs(10)

@everywhere using Distributions, Clustering, NearestNeighbors, Random, Plots, LaTeXStrings, Colors
@everywhere using Measures
@everywhere using Combinatorics, Graphs, LinearAlgebra
@everywhere using ProgressMeter

@everywhere using SparseArrays
# functions for data handling, like computing the biadjacency and pairwise projections of a hypergraphs from it's edges 
@everywhere function get_biadj(hedges,n)
    # n - number of nodes
    # m - number of hyperedges
    # this function returns the biadjacency matrix "biadj" of size m * n
    # its entry biadj[i,j] = 1 indicates node j contained in hyperedge i

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

# TODO convert this to a sparse function 
@everywhere function get_global_cc_hg(biadj::Matrix)
    # this function computes the global clustering coefficient when we view the hypergraph as a bipartite graph
    # it is same as the function robins_alexander_clustering in NetworkX
    # Robins and Alexander defined bipartite clustering coefficient as four times the number of four cycles C_4 divided by the number of three paths L_3 in a bipartite graph
    m = size(biadj, 1) # number of hyperedges
    edges_size = sum(biadj, dims=2) # edges_size[i] indicates the size of hyperedge i
    common_neighbors_mat = biadj * biadj' # common_neighbors_mat[i,j] = common_neighbors_mat[j,i] indicates the intersection size of hyperedges i and j
    diff_neighbors_mat = zeros(Int, m, m) # diff_neighbors_mat[i,j] is the number of nodes contained in hyperedge i but not in hyperedge j
    for i = 1:m
        for j = 1:m
            diff_neighbors_mat[i,j] = edges_size[i] - common_neighbors_mat[i,j] 
        end
    end

    n4 = 0 # number of 4-cycles (butterflies)
    n3 = 0 # number of 3-paths
    for i = 1:m
        for j = i+1:m
            delta = binomial(common_neighbors_mat[i,j], 2) 
            n4 += delta
            n3 += common_neighbors_mat[i,j] * (diff_neighbors_mat[i,j] + diff_neighbors_mat[j,i]) + 4 * delta
        end
    end
    global_cc_hg = 4 * n4 / n3 
    return global_cc_hg
end

@everywhere function get_global_cc_hg(biadj::SparseMatrixCSC)
    # this function computes the global clustering coefficient when we view the hypergraph as a bipartite graph
    # it is same as the function robins_alexander_clustering in NetworkX
    # Robins and Alexander defined bipartite clustering coefficient as four times the number of four cycles C_4 divided by the number of three paths L_3 in a bipartite graph
    m = size(biadj, 1) # number of hyperedges
    edges_size = sum(biadj, dims=2) # edges_size[i] indicates the size of hyperedge i
    common_neighbors_mat = biadj * biadj' # common_neighbors_mat[i,j] = common_neighbors_mat[j,i] indicates the intersection size of hyperedges i and j
    # diff_neighbors_mat = zeros(Int, m, m) # diff_neighbors_mat[i,j] is the number of nodes contained in hyperedge i but not in hyperedge j
    # for i = 1:m
    #     for j = 1:m
    #         diff_neighbors_mat[i,j] = edges_size[i] - common_neighbors_mat[i,j] 
    #     end
    # end
    #build a function that computes the above dynamically without allocating the memory for it since above is a dense matrix
    function _get_diff_neighbors(i,j)
        return edges_size[i] - common_neighbors_mat[i,j]
    end

    n4 = 0 # number of 4-cycles (butterflies)
    n3 = 0 # number of 3-paths
    # for i = 1:m
    #     for j = i+1:m
    #         delta = binomial(common_neighbors_mat[i,j], 2) 
    #         n4 += delta
    #         n3 += common_neighbors_mat[i,j] * (diff_neighbors_mat[i,j] + diff_neighbors_mat[j,i]) + 4 * delta
    #     end
    # end
    # if common_neighbors_mat[i,j] = 0, we can skip the above 
    rows = rowvals(common_neighbors_mat)
    for j = 1:m 
        i_vals = rows[nzrange(common_neighbors_mat,j)]
        for i in i_vals
            if j>i
                delta = binomial(common_neighbors_mat[i,j],2)
                n4 += delta
                n3 += common_neighbors_mat[i,j] * (_get_diff_neighbors(i,j) + _get_diff_neighbors(j,i)) + 4 * delta
            end
        end
    end
    global_cc_hg = 4 * n4 / n3 
    return global_cc_hg
end

@everywhere function get_global_cc(hedges)
    g = SimpleGraph(pairwise_proj(hedges))
    return global_clustering_coefficient(g)
end

# @everywhere function get_global_cc_g(biadj)
#     # this function projects the hypergraph onto a (unweighted) graph and computes the graph global clustering coefficient
#     # the matrix A is the adj of the projected graph
#     A = biadj' * biadj
#     A .= A .!= 0
#     A[diagind(A)] .= 0
#     g = SimpleGraph(A)
#     global_cc_g = global_clustering_coefficient(g) 
#     return global_cc_g
# end

# TODO run the plotting figures 

Random.seed!(1234)

d = 2
n = 10
n_trials = 25
n_alphas = 25
alphas = range(0,2,length=n_alphas)

@everywhere function run_cc_trial(n,d,alphas)
    X = rand(d,n)

    deg_list = zeros(Int, n)
    degreedist = LogNormal(log(3),1)
    for i = 1:n
        deg_list[i] = min(ceil(Int,rand(degreedist)),n-1)
    end

    # make graphs and convert to biadjacency matrices 
    # println("Computing Biadjacency Matrices...")
    graphs = pmap(alpha -> hypergraph_edges(X, deg_list;radfunc=get_func(alpha)), alphas)
    # println("Computing Robins-Alexander CC")
    global_cc_hg = pmap(x->get_global_cc_hg(SparseMatrixCSC{Int}(get_biadj(x,n))),graphs)
    # println("Computing Global CC")
    global_cc_g = pmap(x->get_global_cc(x),graphs)
    return [global_cc_g, global_cc_hg]
end
function get_cc_row_plotting_data(n=10000,d=2,alphas=range(0,2,25),n_trials=10)
    trials = @showprogress map(x->run_cc_trial(n,d,alphas),1:n_trials)
    global_cc_g_mat = reduce(hcat,first.(trials))
    global_cc_hg_mat = reduce(hcat,map(x->x[2],trials))
    return [global_cc_g_mat, global_cc_hg_mat]
end

println("Simulating data...")
n = 10000
alphas = range(0,2,25)
row1 = get_cc_row_plotting_data(n,2,alphas)
row2 = get_cc_row_plotting_data(n,5,alphas)
row3 = get_cc_row_plotting_data(n,10,alphas)

data = [row1,row2,row3]
row1
# PLOT 1 - quantiles
function make_cc_fig(row_data)
    global_cc_g_mat, global_cc_hg_mat = row_data

    plt = Plots.plot(layout=(1,2), margin=6*Plots.mm)
    Plots.plot!(plt, size=(800,300))

    g_q_info = mapslices(x->quantile(x,[0.0,0.5,1.0]),global_cc_g_mat,dims=2)
    q_lower,q_upper = g_q_info[:,2].-g_q_info[:,1], g_q_info[:,3].-g_q_info[:,2]
    Plots.plot!(alphas, g_q_info[:,2], ribbon=(q_lower,q_upper), 
                leg=false, 
                xlabel=L"\alpha", 
                ylabel="Pairwise CC", 
                guidefontsize = 14,
                linewidth=2, 
                marker=:circle, 
                markerstrokewidth=0, 
                ylims=(0,1),
                subplot=1)

    hg_q_info = mapslices(x->quantile(x,[0.0,0.5,1.0]),global_cc_hg_mat,dims=2)
    q_lower,q_upper = hg_q_info[:,2].-hg_q_info[:,1], hg_q_info[:,3].-hg_q_info[:,2]
    Plots.plot!(alphas, hg_q_info[:,2], ribbon=(q_lower,q_upper),
                leg=false,
                xlabel=L"\alpha", 
                ylabel="Bipartite CC", 
                guidefontsize = 14,
                linewidth=2, 
                marker=:circle, 
                markerstrokewidth=0, 
                ylims=(0,0.5),
                subplot=2)

    #touch up margins
    Plots.plot!(plt[2],left_margin=0Measures.mm)
    return plt 
end

f = make_cc_fig(row1)
Plots.plot!(f,plot_title="n=$n d=2",plot_titlefontsize=20)
f = make_cc_fig(row2)
Plots.plot!(f,plot_title="n=$n d=5",plot_titlefontsize=20)
f = make_cc_fig(row3)
Plots.plot!(f,plot_title="n=$n d=10",plot_titlefontsize=20)
# savefig(plt, "/Users/yuzhu/Desktop/add/n1000.pdf")

