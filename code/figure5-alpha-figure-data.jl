"""
script for showing differences in the choice of epsilon_alpha (radius of DBSCAN)

we show:
 our choice based on an approximate scaling argument 
   if alpha <= 1
    return (dist,deg,dim) -> (alpha^(1/dim)*dist)/(deg)^(1/(2*dim))
  else
    return (dist,deg,dim) -> dist/(deg)^(1/(2*dim)) + (alpha-1)*(dist-dist/(deg)^(1/(2*dim)))
  end
 a linear choice: alpha/2*dist
 a choice without correcting for the dimension: replace alpha^1/dim with alpha
"""

using Distributed
# Distributed.addprocs(10)
@everywhere include("spatial-hypergraph.jl")
@everywhere using Distributions, Random
@everywhere using StatsBase
@everywhere using ProgressMeter
@everywhere using MatrixNetworks, Graphs, SparseArrays
@everywhere using SparseArrays

@everywhere function alpha_func_linear(alpha)
    # (dist,deg,dim) -> (alpha^(1/dim)*dist)/(deg)^(1/(2*dim))
    return (dist,deg,dim) -> alpha/2*dist
end

@everywhere function alpha_func_no_dim(alpha)
    if alpha <= 1
        return (dist,deg,dim) -> (alpha*dist)/(deg)^(1/(2*dim))
    else
        return (dist,deg,dim) -> dist/(deg)^(1/(2*dim)) + (alpha-1)*(dist-dist/(deg)^(1/(2*dim)))
    end
end

@everywhere function run_trial(n,d,alphas,alpha_func = alpha->get_func(alpha))
    X = rand(d,n)

    deg_list = zeros(Int, n)
    degreedist = LogNormal(log(3),1)
    for i = 1:n
        deg_list[i] = min(ceil(Int,rand(degreedist)),n-1)
    end

    graphs = pmap(alpha -> hypergraph_edges(X, deg_list;radfunc=alpha_func(alpha)), alphas)
    nedges = pmap(x->length(x), graphs)
    return nedges
end

function _get_plotting_data(n,d,alphas,ntrials,alpha_func = alpha->get_func(alpha))
    return trials = @showprogress map(x->run_trial(n,d,alphas,alpha_func), 1:ntrials)
end

Random.seed!(112358)
alphas = range(0,2,25)
n = 10000
ntrials = 25
save_data = Dict()
for d in [2,5,10]
    for (alpha_func,alpha_func_name) in zip([get_func,alpha_func_linear,alpha_func_no_dim],["alpha_func","linear","no_dim"])
        trials = _get_plotting_data(n,d,alphas,ntrials,alpha_func)
        nedges = reduce(hcat,trials)
        save_data[(n,d,alpha_func_name)] = nedges
    end
end

# save data 
open("data/output/alpha_func_data-n_$n.json", "w") do file
    JSON.print(file, save_data)
end