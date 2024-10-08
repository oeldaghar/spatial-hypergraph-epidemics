using Distributed
using ProgressMeter
using Random
using JSON

@everywhere include("hypergraph-tools.jl")
@everywhere include("spatial-hypergraph.jl")

#parameters
n = 50000
d = 2
ntrials = 25 
alphas = range(0,2,length=25)

random_seed = 1234
Random.seed!(random_seed)

# function designed to interpolate between pure hypergraph at 2
# and pure graph at 0.
#get_func(alpha) = (d,deg) -> (d/2 - d/sqrt(deg))*alpha^2 + (2*d/sqrt(deg) - d/2)*alpha
@everywhere function get_func(alpha) 
  if alpha <= 1
    return (d,deg) -> alpha*d/sqrt(deg)
  else
    return (d,deg) -> d/sqrt(deg) + (alpha-1)*(d-d/sqrt(deg))
  end
end 

function run_trial()
  X = rand(d,n)
  degreedist=LogNormal(log(3),1)
  degs = rand(degreedist,n)
  degs = min.(ceil.(Int,degs),n-1)

  graphs = pmap(alpha -> hypergraph_edges(X;degs=degs,radfunc=get_func(alpha))[1], alphas)
  nedges = map(x->length(x), graphs)
  return nedges
end
trials = @showprogress map(x->run_trial(), 1:25)  

## plotting average number of edge s 

@everywhere function _project_graph_avgdegrees(edges;scalefunc=x->1/x)
  sumdegs = 0.0
  for e in edges
    sumdegs += length(e)*(length(e)-1)*scalefunc(length(e))
  end
  return sumdegs
end
function run_projected_trial(;scalefunc=x->1/x)
  X = rand(d,n)
  degreedist=LogNormal(log(3),1)
  degs = rand(degreedist,n)
  degs = min.(ceil.(Int,degs),n-1)

  graphs = pmap(alpha -> hypergraph_edges(X;degs=degs,radfunc=get_func(alpha))[1], alphas)
  nedges = map(x->_project_graph_avgdegrees(x;scalefunc), graphs)
  return nedges   
end 

## g(m) = m^2
Random.seed!(random_seed)
projected_trials_squared = @showprogress map(x->run_projected_trial(scalefunc=x->1/x^2), 1:25)  
open("data/output/spatial_graph_projected_degree_$(n)_$(d)_squared.json", "w") do io
  output = Dict("n" => n, "d" => d, 
    "alphas" => collect(alphas),  "projected_trials" => projected_trials_squared, 
    "ntrials" => length(projected_trials_squared),"random_seed"=>random_seed)
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
projected_trials = @showprogress map(x->run_projected_trial(), 1:25)  
open("data/output/spatial_graph_projected_degree_$(n)_$(d)_linear.json", "w") do io
  output = Dict("n" => n, "d" => d, 
    "alphas" => collect(alphas),  "projected_trials" => projected_trials, 
    "ntrials" => length(projected_trials),"random_seed"=>random_seed)
  JSON.print(io, output)
end
