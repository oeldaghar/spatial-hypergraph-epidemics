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

@everywhere function get_local_cc(hedges)
    g = SimpleGraph(pairwise_proj(hedges))
    return mean(local_clustering_coefficient(g))
end

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
    local_cc_g = pmap(x->get_local_cc(x),graphs)
    return [global_cc_g, local_cc_g, global_cc_hg]
end
function get_cc_row_plotting_data(n=10000,d=2,alphas=range(0,2,25),n_trials=10)
    trials = @showprogress map(x->run_cc_trial(n,d,alphas),1:n_trials)
    global_cc_g_mat = reduce(hcat,first.(trials))
    local_cc_g_mat = reduce(hcat,map(x->x[2],trials))
    global_cc_hg_mat = reduce(hcat,map(x->x[3],trials))
    return [global_cc_g_mat, local_cc_g_mat, global_cc_hg_mat]
end

println("Simulating data...")
Random.seed!(1234)
n = 1000
alphas = range(0,2,25)
row1 = get_cc_row_plotting_data(n,2,alphas)
row2 = get_cc_row_plotting_data(n,5,alphas)
row3 = get_cc_row_plotting_data(n,10,alphas)

data = [row1,row2,row3]
# PLOT 1 - quantiles
function quantile_plot(data_mat,alphas)
    q_info = mapslices(x->quantile(x,[0.0,0.5,1.0]),data_mat,dims=2)
    q_lower = q_info[:,2].-q_info[:,1]
    q_upper = q_info[:,3].-q_info[:,2]
    f = Plots.plot(alphas, q_info[:,2], ribbon=(q_lower,q_upper), 
                leg=false, 
                # xlabel=L"\alpha", 
                # ylabel="Pairwise CC", 
                guidefontsize = 14,
                linewidth=2, 
                marker=:circle, 
                markerstrokewidth=0, 
                ylims=(0,1),
    )   
    return f 
end


function make_cc_fig(row_data)
    global_cc_g_mat, local_cc_g_mat ,global_cc_hg_mat = row_data

    f1 = quantile_plot(global_cc_g_mat,alphas)
    Plots.plot!(f1,
                xlabel=L"\alpha", 
                ylabel="Global Pairwise CC", 
            )

    f2 = quantile_plot(local_cc_g_mat,alphas)
    Plots.plot!(f2,
                xlabel=L"\alpha", 
                ylabel="Local Pairwise CC", 
                ylims = (0.2,0.85)
            )

    f3 = quantile_plot(global_cc_hg_mat,alphas)
    Plots.plot!(f3,
                xlabel=L"\alpha", 
                ylabel=L"CC_4", 
                ylims=(0,0.35)
            )
    plt = Plots.plot(f1,f2,f3,layout=(1,3), margin=6*Plots.mm,size=(1800,400))
                                        
    #touch up margins
    Plots.plot!(plt[2],left_margin=2Measures.mm,bottom_margin=5Measures.mm)
    Plots.plot!(plt[1],left_margin=10Measures.mm)
    return plt 
end

f1 = make_cc_fig(row1)
Plots.plot!(f1,plot_title="n=$n d=2",plot_titlefontsize=20)
f2 = make_cc_fig(row2)
Plots.plot!(f2,plot_title="n=$n d=5",plot_titlefontsize=20)
f3 = make_cc_fig(row3)
Plots.plot!(f3,plot_title="n=$n d=10",plot_titlefontsize=20)

# put them all together
Plots.plot(f1,f2,f3,layout=(3,1),size=(1800,1200),
        top_margin=-5Plots.mm)


# Tracking the evolution of a single node whose local CC decreases across alpha

@everywhere function track_local_cc(n,d,alphas)
    X = rand(d,n)

    deg_list = zeros(Int, n)
    degreedist = LogNormal(log(3),1)
    for i = 1:n
        deg_list[i] = min(ceil(Int,rand(degreedist)),n-1)
    end

    # make graphs and convert to biadjacency matrices 
    # println("Computing Biadjacency Matrices...")
    graphs = pmap(alpha -> hypergraph_edges(X, deg_list;radfunc=get_func(alpha)), alphas)
    pairwise_graphs = pmap(x->pairwise_proj(x),graphs)
    
    local_cc_g = pmap(x->local_clustering_coefficient(SimpleGraph(x)),pairwise_graphs)
    return X,deg_list,graphs,pairwise_graphs,local_cc_g
end

res = track_local_cc(1000,2,range(0,2,10))

reduce(hcat,res[end])

res[4][end][:,1]