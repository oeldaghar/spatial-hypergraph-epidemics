using Distributions, Clustering, NearestNeighbors
function hypergraph_edges(X;degreedist,radfunc=(dist,deg) -> dist/sqrt(deg))
  T = BallTree(X)
  # form the edges for sparse
  edges = Vector{Int}[]
  for i in axes(X,2)
    #deg = min(ceil(Int,rand(degreedist)),n-1)
    deg = degreedist[i]
    idxs, dists = knn(T, X[:,i], deg+1, true)
    if deg > 1 
      maxdist = maximum(dists) 
      pts = @view X[:,idxs]
      rad = radfunc(maxdist,deg)
      clusters = dbscan(pts, rad).clusters
      for c in clusters
        e = [i]
        for v in c.core_indices
          if idxs[v] != i
            push!(e, idxs[v])
          end
        end
        for v in c.boundary_indices
          if idxs[v] != i 
            push!(e, idxs[v])
          end
        end
        if length(e) > 1
          push!(edges, e)
        end
      end
    else
      # only one vertex! 
      push!(edges, [i,idxs[2]])
    end 
  end 
  return edges, X
end

function spatial_hypergraph_edges(n::Integer,d::Integer;kwargs...)
  X = rand(d,n)
  return hypergraph_edges(X;kwargs...)
end

function get_func(alpha) 
  if alpha <= 1
    return (d,deg) -> alpha*d/sqrt(deg)
  else
    return (d,deg) -> d/sqrt(deg) + (alpha-1)*(d-d/sqrt(deg))
  end
end 

#using Arpack, MatrixNetworks, SparseArrays, LinearMaps
using GenericArpack, MatrixNetworks, SparseArrays

function _hypergraph_components(edges)
  # compute the bipartite graph for the hypergraph adjacency matrix 
  n = maximum(map(x->maximum(x), edges))
  m = length(edges)
  bsrc = Int[]
  bdst = Int[]
  for (i,e) in pairs(edges)
    for v in e
      push!(bsrc, i+n)
      push!(bdst, v)
      push!(bsrc, v)
      push!(bdst, i+n)
    end
  end

  A = sparse(bsrc,bdst,1.0,m+n,m+n)

  scc = scomponents(A)
  return length(scc.sizes) # this is the number of connected components
end 
function _hypergraph_largest_component(edges)
  # compute the bipartite graph for the hypergraph adjacency matrix 
  n = maximum(map(x->maximum(x), edges))
  m = length(edges)
  bsrc = Int[]
  bdst = Int[]
  for (i,e) in pairs(edges)
    for v in e
      push!(bsrc, i+n)
      push!(bdst, v)
      push!(bsrc, v)
      push!(bdst, i+n)
    end
  end

  A = sparse(bsrc,bdst,1.0,m+n,m+n)

  scc = scomponents(A)
  bigcomp = argmax(scc.sizes)

  # now we need to find the edges and nodes that are in the big component. 
  nodes = Int[] 
  bigcc_edges_ids = Set{Int}()
  for (i,c) in pairs(scc.map) 
    if c == bigcomp
      if i <= n
        push!(nodes, i)
      else
        push!(bigcc_edges_ids, i-n)
      end
    end
  end

  # Now label the nodes in the big component
  nodemap = zeros(n) 
  sort!(nodes) # sort the list of nodes. 
  for (i,v) in pairs(nodes)
    nodemap[v] = i
  end

  # Now we rewrite the edges in the large component ...
  newedges = Vector{Int}[]
  for i in bigcc_edges_ids
    e = edges[i]
    newe = map(x->nodemap[x], e)

    @assert all(x->x>0, newe) # make sure all the nodes are in the big component
    push!(newedges, newe)
  end

  return newedges 
end 
function _projected_graph_lam1(edges;scalefunc=x->1/x)
  n = maximum(map(x->maximum(x), edges))

  # create a function to compute A*x where A is the projected graph
  function _matvec!(y,x)
    fill!(y,0.0)
    for e in edges
      edgescale = scalefunc(length(e)-1)
      # we are going to compute the rank-1 update for each edge
      # which is (e*e'*edgescale - I*edgescale) for everything in the edge
      esum = 0.0
      for i in e 
        esum += x[i]*edgescale 
      end 
      for i in e
        y[i] += esum - x[i]*edgescale
      end
    end
    return y 
  end 

  # Amat = zeros(n,n)
  # for i=1:n
  #   x = zeros(n)
  #   x[i] = 1.0
  #   y = zeros(n)
  #   _matvec!(y,x)
  #   Amat[:,i] .= y
  # end
  # display(Amat)

  fop = ArpackSimpleFunctionOp((y,x) -> _matvec!(y,x), n)
  vals, _ = symeigs(fop, 4; which=:LM, maxiter=10000) 
  
  sort!(vals; rev=true, by=abs)
  #@show vals 
  return vals[1]
end

function _project_graph_matrix(edges;scalefunc=x->1/x)
  n = maximum(map(x->maximum(x), edges))
  A = zeros(n,n)
  for e in edges 
    for i in e
      for j in e
        if i != j
          A[i,j] += scalefunc(length(e)-1)
        end
      end
    end
  end
  return A
end

d = 2
n = 500
using Random, LinearAlgebra
Random.seed!(1234)
X = rand(d,n)
degs = min.(ceil.(Int,rand(LogNormal(log(3),1),size(X,2))),size(X,2)-1)
alphas = range(0,2,length=25)
graphs = map(alpha -> hypergraph_edges(X;degreedist=degs,radfunc=get_func(alpha))[1], alphas)

function test_lam1_vs_matrix(edges)
  A = _project_graph_matrix(edges)
  vals = eigvals(A) 
  sort!(vals; rev=true, by=abs)
  lam1 = vals[1] 
  #@show vals[1:4]

  vals2, _ = symeigs(A, 4; which=:LM, maxiter=10000)
  sort!(vals2; rev=true, by=abs)

  #@show vals2[1:4]

  lam1_matvec = _projected_graph_lam1(edges)
  return lam1, lam1_matvec
end 

test_lam1_vs_matrix(graphs[2])

## 
using Test 
for (i,g) in pairs(graphs)
  vals = test_lam1_vs_matrix(g)
  println("$i: ", vals)
  @test vals[1] ≈ vals[2]
end



## Now let's plot the average number of edges...
using Plots
function plot_data(trials)
  nedges = reduce(hcat, trials) # now we have alpha as row index and trial as column index

  p = Plots.plot(xlabel="α", ylabel="Number of edges")

  # show the full data as a ribbon
  minl = minimum.(eachrow(nedges))
  maxl = maximum.(eachrow(nedges))
  mid = (minl .+ maxl)/2
  upper = maxl .- mid
  lower = mid .- minl
  @show minl, maxl
  Plots.plot!(p, alphas, (minl .+ maxl)/2, linewidth=0, ribbon=(upper, lower), label="", c=Gray(0.8))

  # I want to take the quantiles over the columns... 
  quantiles = [0.1,0.5,0.9]
  linewidths = [1.5, 3, 1.5]
  colors = [ :grey, :black, :grey,]

  for (q, lw, c) in zip(quantiles, linewidths, colors)
    nedges_q = quantile.(eachrow(nedges), q)
    Plots.plot!(p, alphas, nedges_q, label="", linewidth=lw, color=c)
  end
  
  p
end

## Now I want to project each graph
using Random, JSON, ProgressMeter
n = 50000
d = 2 
random_seed = 1234 
ntrials = 25 
degreedist = LogNormal(log(3),1)
Random.seed!(random_seed)
seeds = rand(UInt, ntrials )

function run_projected_trial(seed;scalefunc=x->1/x)
  Random.seed!(seed)
  X = rand(d,n)
  degs = rand(degreedist, n)
  degs = min.(ceil.(Int,degs),n-1)
  graphs = @showprogress offset=2 desc="Make graphs ..." map(alpha -> hypergraph_edges(X;degreedist=degs, radfunc=get_func(alpha))[1], alphas)
  nedges = @showprogress offset=2 desc="Eigenvalues ..." map(x->_projected_graph_lam1(_hypergraph_largest_component(x);scalefunc), graphs)
  return nedges   
end 

##
using ProgressMeter
projected_trials = @showprogress offset=1 desc="Run trials ..." map(x->run_projected_trial(x), seeds)  

##
open("spatial_graph_projected_lambda1_$(n)_$(d)_linear.json", "w") do io
  output = Dict("n" => n, "d" => d, 
    "degreedist_mu" => degreedist.μ, "degreedist_sigma" => degreedist.σ,
    "alphas" => collect(alphas),  "projected_trials" => projected_trials, 
    "ntrials" => length(projected_trials), "random_seed" => random_seed, 
    "seeds" => seeds)
  JSON.print(io, output)
end

##
using LaTeXStrings
p = plot_data(projected_trials)
plot!(p, ylabel="Weighted projected λ₁", xlabel="α",size=(300,300))
p

##
savefig(p, "spatial_graph_projected_lambda1_$(n)_$(d)_linear.pdf")




##
projected_trials = @showprogress offset=1 desc="Run trials ..." map(x->run_projected_trial(x;scalefunc=x->1/x^2), seeds)  

##
open("spatial_graph_projected_lambda1_$(n)_$(d)_squared.json", "w") do io
  output = Dict("n" => n, "d" => d, 
    "degreedist_mu" => degreedist.μ, "degreedist_sigma" => degreedist.σ,
    "alphas" => collect(alphas),  "projected_trials" => projected_trials, 
    "ntrials" => length(projected_trials), "random_seed" => random_seed, 
    "seeds" => seeds)
  JSON.print(io, output)
end

##
using LaTeXStrings
p = plot_data(projected_trials)
plot!(p, ylabel="Weighted projected λ₁", xlabel="α",size=(300,300))
p

##
savefig(p, "spatial_graph_projected_lambda1_$(n)_$(d)_squared.pdf")

##
projected_trials = @showprogress offset=1 desc="Run trials ..." map(x->run_projected_trial(x;scalefunc=x->1/sqrt(x)), seeds)  

##
open("spatial_graph_projected_lambda1_$(n)_$(d)_sqrt.json", "w") do io
  output = Dict("n" => n, "d" => d, 
    "degreedist_mu" => degreedist.μ, "degreedist_sigma" => degreedist.σ,
    "alphas" => collect(alphas),  "projected_trials" => projected_trials, 
    "ntrials" => length(projected_trials), "random_seed" => random_seed, 
    "seeds" => seeds)
  JSON.print(io, output)
end

##
using LaTeXStrings
p = plot_data(projected_trials)
plot!(p, ylabel="Weighted projected λ₁", xlabel="α",size=(300,300))
p

##
savefig(p, "spatial_graph_projected_lambda1_$(n)_$(d)_sqrt.pdf")