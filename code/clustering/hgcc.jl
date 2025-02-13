using Distributed 
# addprocs(10)

@everywhere include("../spatial-hypergraph.jl")
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
    fill!(A.nzval,1.0)
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

println("Simulating data...")
Random.seed!(1234)
n = 1000
alphas = range(0,2,25)
row1 = get_cc_row_plotting_data(n,2,alphas)
row2 = get_cc_row_plotting_data(n,5,alphas)
row3 = get_cc_row_plotting_data(n,10,alphas)


row1[3]

data = [row1,row2,row3]

f1 = make_cc_fig(row1)
Plots.plot!(f1,plot_title="n=$n d=2",plot_titlefontsize=20)
f2 = make_cc_fig(row2)
Plots.plot!(f2,plot_title="n=$n d=5",plot_titlefontsize=20)
f3 = make_cc_fig(row3)
Plots.plot!(f3,plot_title="n=$n d=10",plot_titlefontsize=20)

# put them all together
plt = Plots.plot(f1,f2,f3,layout=(3,1),size=(1400,1200),
        top_margin=-5Plots.mm)

Plots.savefig(plt,"data/output/figures/final/hypergraph-cc.pdf")

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


Random.seed!(278) #small intermediate decrease
alphas = range(0,2,10)
X,deg_list,graphs,pairwise_graphs,local_cc_g = track_local_cc(200,2,alphas);

dvecs = reduce(hcat,map(x->vec(sum(x,dims=1)),pairwise_graphs));
local_cc_g = reduce(hcat,local_cc_g);
# all elements that experience a decrease 
inds = findall((local_cc_g[:,2].-local_cc_g[:,1]).<0)
val,ind = findmin(dvecs[inds,end])

node_id = inds[ind]
dvecs[node_id,end]
# find first decrease for this node 
first_decrease = findfirst([local_cc_g[node_id,i+1]-local_cc_g[node_id,i]<0 for i=1:lastindex(local_cc_g,2)-1])+1

local_cc_g[node_id,first_decrease-1]
local_cc_g[node_id,first_decrease]

function get_ego_info(node,hedges,X)
    # find all hyperedges containing our node 
    hedges_ind = findall([node in x for x in hedges])
    ego_hedges = hedges[hedges_ind]
    ego_neighbors = reduce(vcat,ego_hedges)
    ego_neighbors = sort(unique(ego_neighbors))
    
    # induced_edges = []
    # for edge in hedges
    #     if any([v in edge for v in ego_neighbors])
    #         push!(induced_edges,edge)
    #     end
    # end
    # induced_nodes = sort(unique(reduce(vcat,induced_edges)))
    # get data coordinates for ego neighbors 
    return ego_hedges,ego_neighbors,X[:,ego_neighbors] #, induced_edges, induced_nodes, X[:,induced_nodes]
end

function get_induced_info(nodes,hedges,X)
    # find all hyperedges containing our node 
    induced_edges = []
    for edge in hedges 
        if any([node in edge for node in nodes])
            push!(induced_edges,edge)
        end
    end
    induced_nodes = sort(unique(reduce(vcat,induced_edges)))
    return induced_edges, induced_nodes
end


# get_ego_info(node_id,graphs)
# get_induced_info(old_neighbors,graphs[1],X)[2]

# pairwise projections 
using LazySets
using MatrixNetworks

function plot_hypergraph(X,hedges,nodes,projected=false)
    plt = Plots.plot()
    for edge in hedges
        if length(edge) > 2 && projected==false
            local b = [Ball2(X[:, v], 0.0005) for v in edge]
            local c = ConvexHullArray(b)
            Plots.plot!(plt, c, 1e-3, alpha=0.05, c=:blue)
        else
            Plots.plot!(plt, X[1,edge], X[2,edge], c=:blue)
        end
    end
    xs = [X[1,v] for v in nodes]
    ys = [X[2,v] for v in nodes]
    Plots.scatter!(plt,xs,ys,leg=false,markerstrokewidth=0,markersize=12,markercolor=:blue)
    return plt
end


figs = []
for alpha_ind = 1:lastindex(graphs)
    # hypergraph vis 1
    curr_hedge,curr_neighbors,curr_X = get_ego_info(node_id,graphs[alpha_ind],X)
    induced_edges,induced_nodes = get_induced_info(curr_neighbors,graphs[alpha_ind],X)

    f = plot_hypergraph(X, induced_edges, curr_neighbors)
    nwedges = binomial(length(curr_neighbors)-1,2)
    ntris = length(collect(MatrixNetworks.triangles(pairwise_graphs[alpha_ind],node_id)))
    ratio = round(ntris/nwedges,digits=2)
    alpha_val = round(alphas[alpha_ind],digits=2)
    Plots.plot!(f,annotation=(0.46,0.12,"alpha: $(alpha_val)\nWedges: $nwedges\nTriangles: $ntris\nRatio: $ratio"))
    Plots.scatter!(f,[X[1,node_id]],[X[2,node_id]],markercolor=:red,markersize=12,
            ylims = (0.05,0.16),
            xlims = (0.43,0.47)
    )
    push!(figs,deepcopy(f))
end
Plots.plot(figs...,layout=(4,3),size=(600*3,400*4),link=:all)



#### non-induced edges
figs = []
for alpha_ind = 1:lastindex(graphs)
    # hypergraph vis 1
    curr_hedge,curr_neighbors,curr_X = get_ego_info(node_id,graphs[alpha_ind],X)
    # induced_edges,induced_nodes = get_induced_info(curr_neighbors,graphs[alpha_ind],X)

    f = plot_hypergraph(X, curr_hedge, curr_neighbors)
    nwedges = binomial(length(curr_neighbors)-1,2)
    ntris = length(collect(MatrixNetworks.triangles(pairwise_graphs[alpha_ind],node_id)))
    ratio = round(ntris/nwedges,digits=2)
    alpha_val = round(alphas[alpha_ind],digits=2)
    Plots.plot!(f,annotation=(0.46,0.12,"alpha: $(alpha_val)\nWedges: $nwedges\nTriangles: $ntris\nRatio: $ratio"))
    Plots.scatter!(f,[X[1,node_id]],[X[2,node_id]],markercolor=:red,markersize=12,
            ylims = (0.05,0.16),
            xlims = (0.43,0.47)
    )
    push!(figs,deepcopy(f))
end
Plots.plot(figs...,layout=(4,3),size=(600*3,400*4),link=:all)





function commonNeighbors(A::SparseMatrixCSC)
    @assert issymmetric(A)
    P = deepcopy(A)
    fill!(nonzeros(P),0)
    rowval = rowvals(A)
    neighs = zeros(Bool, size(A,1))
    for j=1:size(A,1)
        # index neighbors
        for nzj in nzrange(A, j)
            neighs[rowval[nzj]] = true
        end
        # for each edge out of this node...
        for nzj in nzrange(A, j)
            i = rowval[nzj]
            score = 0
            for nzi in nzrange(A, i)
                w = rowval[nzi]
                if neighs[w] == true
                    # for each neighbor of i (w) that is
                    # also a neighbor of j neighs[w] = true
                    # increase the score
                    score += 1
                end
            end
            nonzeros(P)[nzj] = score
        end
        # reset indicies
        for nzj in nzrange(A, j)
            neighs[rowval[nzj]] = false
        end
    end
    return P
end


function multigraph_triangles(A::SparseMatrixCSC)
    @assert(isequal(size(A)...))
    rowval = rowvals(A)
    d = vec(sum(A;dims=1))
    neighs = zeros(Float64,lastindex(A,1))
    wedge_counts = zeros(Int,lastindex(A,1))
    tri_counts = zeros(Int,lastindex(A,1))
    for v = 1:lastindex(A,1)
        # cache neighbors and their weights 
        for u in rowval[nzrange(A,v)]
            neighs[u] = A[u,v]
        end
        # for each neighbor, iterate over neighbors (1-step BFS) and track triangles
        for u in rowval[nzrange(A,v)]
            edge_weight = A[u,v]
            for w in rowval[nzrange(A,u)]
                if w!=v
                    possible_tris = edge_weight*A[w,u]
                    if neighs[w]>0
                        tri_counts[u] += possible_tris
                    end
                    wedge_counts[u] += possible_tris
                end
            end
        end
        # reset neighbors 
        for u in rowval[nzrange(A,v)]
            neighs[u] = 0.0
        end
    end
    return tri_counts,wedge_counts
end

A = pairwise_proj(graphs[1])

dropzeros!(A)
unique(A.nzval)

G = SimpleGraph(A)
local_clustering_coefficient(G)
tri,weds = multigraph_triangles(A)

norm(tri./weds .- local_clustering_coefficient(G))
norm(Graphs.triangles(G)*2-tri)

Graphs.triangles(G)./local_clustering_coefficient(G)
weds


A = sparse([0 0 0 1; 0 0 0 1; 0 0 0 1;1 1 1 0])
multigraph_triangles(A)
Graphs.triangles(SimpleGraph(A))
