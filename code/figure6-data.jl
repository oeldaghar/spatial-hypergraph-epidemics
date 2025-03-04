using Distributed
using JSON
# nworkers()
# Distributed.addprocs(10)
@everywhere include("spatial-hypergraph.jl")
@everywhere using Distributions, Random
@everywhere using StatsBase
@everywhere using ProgressMeter
@everywhere using MatrixNetworks, Graphs, SparseArrays
@everywhere using SparseArrays


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

@everywhere function run_trial(n,d,alphas)
    X = rand(d,n)
    
    deg_list = zeros(Int, n)
    degreedist = LogNormal(log(3),1)
    for i = 1:n
        deg_list[i] = min(ceil(Int,rand(degreedist)),n-1)
    end

    graphs = pmap(alpha -> hypergraph_edges(X, deg_list;radfunc=get_func(alpha)), alphas)
    
    nedges = pmap(x->length(x), graphs)
    ntris = pmap(x->sum(1 for i in MatrixNetworks.triangles(pairwise_proj(x))), graphs)
    
    edge_cnts = pmap(x->get_edge_counts(x),graphs)
    
    return nedges,ntris,edge_cnts
end

function get_row_plotting_data(n=10000,d=2,alphas=range(0,2,25),ntrials=25)

    trials = @showprogress map(x->run_trial(n,d,alphas), 1:ntrials)  

    nedges = reduce(hcat, first.(trials)) # now we have alpha as row index and trial as column index
    ntris = reduce(hcat, map(x->x[2],trials))
    # handle edge_cnts. comes in as a vector for each alpha and trial
    # matrix of alphas vs trial. vector in each component
    edge_cnts = reduce(hcat,map(x->x[3],trials))
    # aggreagte across trials #ncols 
    edge_cnts = sum_columns_of_matrix(edge_cnts)./ntrials
    # edge_cnts = edge_cnts./=sum(edge_cnts,dims=1)

    return [edge_cnts, nedges, ntris]
end 


alphas = range(0,2,25)
n = 10000
ntrials = 25
Random.seed!(457)
println("Simulating data...")
row1 = get_row_plotting_data(n,2,alphas,ntrials)
row2 = get_row_plotting_data(n,5,alphas,ntrials)
row3 = get_row_plotting_data(n,10,alphas,ntrials)
data = [row1,row2,row3]


# make save data 
save_data = Dict()
for (row,key) in zip(data,[(n,2),(n,5),(n,10)])
    for (mat,mat_name) in zip(row,["edge_cnts", "nedges", "ntris"])
        save_data[tuple(key...,mat_name)] = mat
    end
end
open("data/output/hypergraph_stats-n_$n.json", "w") do file
    JSON.print(file, save_data)
end
