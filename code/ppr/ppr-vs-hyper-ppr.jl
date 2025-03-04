using Distributed
# Distributed.addprocs(10)
@everywhere include("../spatial-hypergraph.jl")
@everywhere using Distributions, Random
@everywhere using StatsBase
@everywhere using ProgressMeter
@everywhere using MatrixNetworks, Graphs, SparseArrays
@everywhere using SparseArrays

using JSON

# uses Meng's code from https://github.com/MengLiuPurdue/LHQD/tree/main
include("LHQD/local-hyper.jl")
include("LHQD/common.jl")
### EXAMPLE FROM MENG DOCS
# H should be biadjacency matrix of the form E X V (SparseCSC)
# G = LH.graph(H,1.0) # a hypergraph object with delta=1.0 in its cut function
# q = 2.0
# L = LH.loss_type(q) # the loss-type, this is a 2-norm)
# kappa = 0.01 # value of kappa (sparsity regularization)
# gamma = 0.1 # value of gamma (regularization on seed) 
# rho = 0.5 # value of rho (KKT apprx-val)
# S = [1] # seed set
# x,r,iter = LH.lh_diffusion(G,S,gamma,kappa,rho,L,max_iters=10000)
# cond,cluster = hyper_sweepcut(G.H,x,G.deg,G.delta,0.0,G.order)
### END EXAMPLE

# hedges to biadjacency E X V as SparseCSC
function hedges_to_biadj(hedges,n)
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

# hypergraph ppr
function hypergraph_ppr(H,S)
    # H = hedges_to_biadj(hedges,n) # biadjacency rep
    G = LH.graph(H,1.0) # a hypergraph object with delta=1.0 in its cut function
    x,r,iter = LH.lh_diffusion(G,S,gamma,kappa,rho,L,max_iters=10000)
    # cond,cluster = hyper_sweepcut(G.H,x,G.deg,G.delta,0.0,G.order)
    return x 
end


n,d = 500,2
alphas = range(0,2,3)

# hypergraph ppr parameters 
q = 2.0
L = LH.loss_type(q) # the loss-type, this is a 2-norm)
kappa = 0.0001 # value of kappa (sparsity regularization)
gamma = 1.0 # value of gamma (regularization on seed) 
rho = 0.5 # value of rho (KKT apprx-val)
# S = [rand(1:n)] # seed set

Random.seed!(11)
X = rand(d,n)

# get top right corner node for visualizations
_,S = findmin(vec(mapslices(xy->(xy[1]-1)^2+(xy[2]-1)^2,X,dims=1)))
S = [S]

deg_list = zeros(Int, n)
degreedist = LogNormal(log(3),1)
for i = 1:n
    deg_list[i] = min(ceil(Int,rand(degreedist)),n-1)
end

# make graphs and convert to biadjacency matrices 
graphs = pmap(alpha -> hypergraph_edges(X, deg_list;radfunc=get_func(alpha)), alphas)

ppr_soln = map(x->hypergraph_ppr(hedges_to_biadj(x,n),S),graphs)
ppr_soln = reduce(hcat,ppr_soln)


# save edges, spatial info, and higher order ppr solutions 
save_data = Dict()
save_data["X"] = X
save_data["graphs"] = graphs
save_data["ppr_soln"] = ppr_soln 
save_data["alphas"] = collect(alphas)

# save data 
open("data/output/hyper_ppr_data.json", "w") do io
    JSON.print(io, save_data)
end