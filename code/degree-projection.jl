using Distributed
using JSON

# addprocs(10)
@everywhere using ProgressMeter
@everywhere using Random

@everywhere include("hypergraph-tools.jl")
@everywhere include("spatial-hypergraph.jl")

#parameters
@everywhere n = 50000
@everywhere d = 2
@everywhere alphas = range(0,2,length=25)
ntrials = 25 

random_seed = 1234
Random.seed!(random_seed)

@everywhere function run_trial()
  X = rand(d,n)
  degreedist=LogNormal(log(3),1)
  degs = rand(degreedist,n)
  degs = min.(ceil.(Int,degs),n-1)

  graphs = map(alpha -> hypergraph_edges(X;degs=degs,radfunc=get_func(alpha))[1], alphas)
  nedges = map(x->length(x), graphs)
  return nedges
end
trials = @showprogress pmap(x->run_trial(), 1:25)  

## plotting average number of edge s 

@everywhere function _project_graph_avgdegrees(edges;scalefunc=x->1/x)
  sumdegs = 0.0
  for e in edges
    sumdegs += length(e)*(length(e)-1)*scalefunc(length(e))
  end
  return sumdegs
end
@everywhere function run_projected_trial(;scalefunc=x->1/x)
  X = rand(d,n)
  degreedist=LogNormal(log(3),1)
  degs = rand(degreedist,n)
  degs = min.(ceil.(Int,degs),n-1)

  graphs = pmap(alpha -> hypergraph_edges(X;degs=degs,radfunc=get_func(alpha))[1], alphas)
  nedges = pmap(x->_project_graph_avgdegrees(x;scalefunc), graphs)
  return nedges   
end 

## g(m) = 1
Random.seed!(random_seed)
projected_trials_pairwise = @showprogress map(x->run_projected_trial(scalefunc=x->1), 1:25)  
open("data/output/spatial_graph_projected_degree_$(n)_$(d)_pairwise.json", "w") do io
  output = Dict("n" => n, "d" => d, 
    "alphas" => collect(alphas),  "projected_trials" => projected_trials_pairwise, 
    "ntrials" => length(projected_trials_pairwise),"random_seed"=>random_seed)
  JSON.print(io, output)
end

## g(m) = sqrt(m)
Random.seed!(random_seed)
projected_trials_sqrt = @showprogress map(x->run_projected_trial(scalefunc=x->1/sqrt(x)), 1:25)  
open("data/output/spatial_graph_projected_degree_$(n)_$(d)_sqrt.json", "w") do io
  output = Dict("n" => n, "d" => d, 
    "alphas" => collect(alphas),  "projected_trials" => projected_trials_sqrt, 
    "ntrials" => length(projected_trials_sqrt),"random_seed"=>random_seed)
  JSON.print(io, output)
end

## g(m)=m
Random.seed!(random_seed)
projected_trials = @showprogress map(x->run_projected_trial(scalefunc=x->1/x), 1:25)  
open("data/output/spatial_graph_projected_degree_$(n)_$(d)_linear.json", "w") do io
  output = Dict("n" => n, "d" => d, 
    "alphas" => collect(alphas),  "projected_trials" => projected_trials, 
    "ntrials" => length(projected_trials),"random_seed"=>random_seed)
  JSON.print(io, output)
end
