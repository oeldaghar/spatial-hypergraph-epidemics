using Distributed 
# addprocs(25)

@everywhere include("../spatial-hypergraph.jl")
@everywhere using Distributions, Clustering, NearestNeighbors, Random
@everywhere using Combinatorics, Graphs, LinearAlgebra
@everywhere using ProgressMeter
@everywhere using MatrixNetworks

# pairwise projections 
using JSON

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

@everywhere function pairwise_proj(hedges;multigraph=false)
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
    A = sparse(ei,ej,ones(Float64,lastindex(ei)),nnodes,nnodes)
    for i=1:nnodes
        A[i,i]=0
    end
    dropzeros!(A)
    if !multigraph # retain edge weights as multigraph representation
        fill!(A.nzval,1.0)
    end
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

@everywhere function multigraph_clustering(A::SparseMatrixCSC)
    # Ensure input matrix is square
    if size(A, 1) != size(A, 2)
        throw(ArgumentError("Input matrix is not square. Got size: $(size(A))"))
    end

    rowval = rowvals(A)
    n = lastindex(A, 1)

    last_seen = zeros(Int, n)       # Counter-based neighbor tracking
    current_iteration = 0          # Iteration counter

    wedge_counts = zeros(Int, n)
    tri_counts = zeros(Int, n)

    neighbors = [rowval[nzrange(A, v)] for v in 1:lastindex(A, 1)]

    @inbounds for v = 1:n
        current_iteration += 1

        # Cache neighbors using counter-based tracking
        @inbounds for u in neighbors[v]
            last_seen[u] = current_iteration
        end

        # Iterate over neighbors and track wedges
        @inbounds for u in neighbors[v]
            @inbounds for w in neighbors[u]
                if w != v
                    possible_wedges = A[u, v] * A[w, u]
                    wedge_counts[u] += possible_wedges
                    if last_seen[w] == current_iteration
                        tri_counts[u] += clamp(A[w,v],0,possible_wedges)
                    end
                end
            end
        end
    end
    # each wedge gets counted twice (one from either side of the center)
    return div.(tri_counts, 2), div.(wedge_counts, 2)
end

@everywhere function multigraph_clustering_coefficient(A::SparseMatrixCSC)
    t,w = multigraph_clustering(A)
    results = zeros(Float64,lastindex(t))
    for ind = 1:lastindex(t)
        if w[ind]>0
            results[ind] = t[ind]/w[ind]
        end
    end
    return results 
end

@everywhere function get_multigraph_cc(hedges)
    A = pairwise_proj(hedges,multigraph=true)
    cc = multigraph_clustering_coefficient(A)
    return mean(cc)
end 

# @everywhere function get_weighted_multigraph_cc(hedges)
#     # computes a weighted clustering coefficient. 
#     # see https://journals.aps.org/pre/pdf/10.1103/PhysRevE.71.065103?casa_token=NGVROFqvWlUAAAAA%3AncYoO3BIiC7LfNFF1IG4HUW1DLMEk0GF2exkNVq59E_BbvL4U1GWxUCwB1gtgZ9DLlc53P22FX-cb47W
#     A = pairwise_proj(hedges,multigraph=true)
#     A.nzval ./= maximum(A.nzval)
#     return mean(MatrixNetworks.clustercoeffs(A))
# end 

## TODO - might be able to improve runtime (map for each cc trial but pmap across trials)

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
    # mulitgraph clustering 
    local_cc_mg = pmap(x->get_multigraph_cc(x),graphs)
    
    return [global_cc_g, local_cc_g, global_cc_hg, local_cc_mg]
end

function get_cc_row_plotting_data(n=10000,d=2,alphas=range(0,2,25),n_trials=25)
    trials = @showprogress map(x->run_cc_trial(n,d,alphas),1:n_trials)
    global_cc_g_mat = reduce(hcat,first.(trials))
    local_cc_g_mat = reduce(hcat,map(x->x[2],trials))
    global_cc_hg_mat = reduce(hcat,map(x->x[3],trials))
    local_cc_mg_mat = reduce(hcat,map(x->x[4],trials))
    return [global_cc_g_mat, local_cc_g_mat, global_cc_hg_mat, local_cc_mg_mat]
end


println("Simulating data...")
Random.seed!(1234)
n = 10000
alphas = range(0,2,25)
row1 = get_cc_row_plotting_data(n,2,alphas)
row2 = get_cc_row_plotting_data(n,5,alphas)
row3 = get_cc_row_plotting_data(n,10,alphas)

data = [row1,row2,row3]

# write to file 
save_data = Dict()
for (row,key) in zip(data,[(n,2),(n,5),(n,10)])
    for (mat_data,mat_name) in zip(row,["global_cc_g_mat", "local_cc_g_mat", "global_cc_hg_mat", "local_cc_mg_mat"])
        save_data[tuple(key...,mat_name)] = mat_data
    end
end
open("data/output/hgcc_data-n_$n.json", "w") do file
    JSON.print(file, save_data)
end


# # star graph 
# A = sparse([0 0 0 1; 0 0 0 1; 0 0 0 1;1 1 1 0])
# @time t,w = multigraph_clustering(A);
# norm(w.-binomial.(vec(sum(A;dims=1)),2))==0
# norm(t .- Graphs.triangles(SimpleGraph(A)))==0

# # random graphs 
# A = sparse(Graphs.erdos_renyi(1000,0.25))
# @time t,w = multigraph_clustering_optimized(A);
# norm(w.-binomial.(vec(sum(A;dims=1)),2))==0
# norm(t .- Graphs.triangles(SimpleGraph(A)))==0

# # specific examples with weights 
# A = sparse([0 2 0;2 0 2;0 2 0])
# t,w = multigraph_clustering_optimized(A)
# norm(t)==0
# all(w .== [0,4,0])

# # specific weighted examples 
# A = sparse([0 2 1;2 0 2;1 2 0])
# t,w = multigraph_clustering_optimized(A)
# norm(t.-[2,1,2])==0
# norm(w .- [2,4,2])==0

# A = sparse([0 2 2;2 0 2;2 2 0])
# t,w = multigraph_clustering_optimized(A)
# norm(t.-[2,2,2])==0
# norm(w .- [4,4,4])==0

# A = sparse([0 2 4;2 0 2;4 2 0])
# t,w = multigraph_clustering_optimized(A)
# norm(t.-[2,4,2])==0
# norm(w .- [8,4,8])==0
# graphs = pmap(alpha -> hypergraph_edges(X, deg_list;radfunc=get_func(alpha)), range(0,2,15))
# get_weighted_multigraph_cc.(graphs)

# function testing_weighted_cc(n,d,ntrials=10)
#     violation = 0
#     @showprogress for trial = 1:ntrials
#         alphas = range(0,2,25)
#         X = rand(d,n)
#         deg_list = zeros(Int, n)
#         degreedist = LogNormal(log(3),1)
#         for i = 1:n
#             deg_list[i] = min(ceil(Int,rand(degreedist)),n-1)
#         end
#         graphs = pmap(alpha -> hypergraph_edges(X, deg_list;radfunc=get_func(alpha)), alphas)
#         weighted_cc = pmap(x->get_weighted_multigraph_cc(x),graphs)
#         if any(isnan.(weighted_cc))
#             println("GOT NAN...")
#             ind = findfirst(isnan.(weighted_cc))
#             return graphs[ind]
#         end
#     end
#     return violation
# end

# res = testing_weighted_cc(500,10,100)
# tmp = deepcopy(res)
# A = pairwise_proj(tmp,multigraph=true)
# A.nzval./=maximum(A.nzval)
# sum(MatrixNetworks.clustercoeffs(A))/lastindex(A,1)

# # write A to file 
# original = MatrixNetworks.clustercoeffs(A)
# for i = 1:100
#     new = MatrixNetworks.clustercoeffs(A)
#     err = norm(new.-original)
#     if err>1e-5
#         println("ERROR, norm difference of $(err)")
#     end
# end

# # writeSMAT 
# function writeSMAT(A::SparseMatrixCSC,filename::String)
#     Is,Js,Vs = findnz(A)
#     nedges = length(Is)
#     open(filename,"w") do file
#         write(file,"$(size(A,1)) $(size(A,2)) $nedges\n")
#         for i=1:nedges
#             #smat assumes node labels start at 0
#             write(file,"$(Is[i]-1) $(Js[i]-1) $(Vs[i])\n") 
#         end
#     end
# end

# writeSMAT(A,"weighted_cc_example.smat")