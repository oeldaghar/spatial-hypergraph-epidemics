##
# generate a random hypergraph

using Random


function _draw_with_duplicate_check(e, rng, elements, duplicates)
  while true
    val = rand(rng, elements) # draw a random elements
    if duplicates == true 
      return val 
    else 
      duplicate = false 
      for curval in e
        if val == curval # found a duplicate 
          duplicate = true 
          break 
        end
      end
      if duplicate
        continue # continue the outer loop
      else
        return val
      end 
    end
  end 
end 
function random_hyperedge!(e, elements, max_size, p; rng, duplicates)
  push!(e, rand(rng, elements))
  push!(e, _draw_with_duplicate_check(e, rng, elements, duplicates))
  for i in 3:max_size
    rv = rand(rng) 
    if rv <= p 
      # then the hyperedge keeps growing... 
      push!(e, _draw_with_duplicate_check(e, rng, elements, duplicates))
    else
      # we stop the hyperedge... 
      break
    end 
  end 
  return e
end 
"""
    random_hyperedge(elements, max_size, p; [rng=GLBOAL_RNG, duplicates=false])

Generate a random hyperedge where the hyperedge keeps growing with probability p. 
(Seeting p=0 will generate a random edge.) The hyperedge will grow up to max_size. 
Elements of the hyperedge are chosen by `rand(rng, elements)`.

We do not allow duplicate elements in the hyperedge if `duplicates` is false. 
**This may cause problems for sets of small elements and cause the algorithm to run forever.**
""" 
function random_hyperedge(elements, max_size::Int, p::Real; rng = Random.GLOBAL_RNG, duplicates::Bool=false)
  e = Vector{eltype(elements)}()
  return random_hyperedge!(e, elements, max_size, p; rng, duplicates)
end 

"""
    random_hypergraph(n, nedges, hep; [rng = Random.GLOBAL_RNG, duplicates=false, max_size])
    random_hypergraph(elements, nedges, hep; [rng = Random.GLOBAL_RNG, duplicates=false, max_size])

Generate a random hypergraph where elements are chosen from 1:n. 
- the default value of max_size is based on twice the expected length of the hyperedge continuation probability.
  So if this is 0.4 (i.e. short hyperedges), then max_size = 2*ceil(Int, 0.4/0.6 ) = 4. 
- duplicates refers to duplicate values within a hyperedge. This function may produce 
  duplicate hyperedges and does not check for uniqueness. If you need that, call 
  `unique!` on the final list of hyperedges. If you need unique sorted hyperedges, then
  use `unique!(sort!.(random_hypergraph()))`
""" 
function random_hypergraph(elements, nedges::Int, hep::Real; rng = Random.GLOBAL_RNG, duplicates::Bool=false, 
  max_size = 2+2*ceil(Int, hep/(1-hep)))

  if duplicates == false && max_size >= length(elements) 
    @warn "This may run forever because max_size is longer than the total list of elements but we don't allow duplicates"
  end 
  edges = Vector{Vector{Int}}() 
  for i=1:nedges
    push!(edges, random_hyperedge(elements, max_size, hep))
  end 
  return edges
end 

#random_hypergraph(n::Int, nedges::Int, hep::Real; kwargs...) = random_hypergraph(1:n, nedges, hep; kwargs... )

#random_hypergraph(10, 50, 0.7)

"""
  planted_partition(n, m, nedges1, nedges2, hep1, hep2; [rng=Random.GLOBAL_RNG, kwargs...])

Generate a planted partition model with n vertices with a planted partition on the first
m of them (region: 1:m).

We generate nedges1 in the full region and nedges2 in the planted region. 
The edges in the full region have hyperedge length probability hep1 and the 
edges in the planted region have length probablitiy hep2 

See random_hypergraph for other kwargs... 
""" 
function planted_partition(n, m, nedges1, nedges2, hep1, hep2; kwargs...)
  m < n || throw(ArgumentError("m must be strictly less than n"))
  edges1 = random_hypergraph(1:n, nedges1, hep1; kwargs...)
  edges2 = random_hypergraph(1:m, nedges2, hep2; kwargs...)
  edges = append!(edges1, edges2) 
  return unique!(sort!.(edges))
end 
# planted_partition(100, 10, 1000, 150, 0.5, 0.5 )

function hypergraph_SBM(n, ncluster, nedges1, nedges2, hep1, hep2; kwargs...) 
    ncluster <= n || throw(ArgumentError("ncluster must be less than n"))
    vlabels = rand(1:ncluster, n)
    sort!(vlabels)
    elabels = rand(1:ncluster, nedges2)
    clusters = Vector{Vector{Int}}()
    edgecnt = zeros(Int64, ncluster)
    edges = random_hypergraph(1:n, nedges1, hep1; kwargs...)
    for i = 1:ncluster
      C = findall(x->vlabels[x]==i, 1:n)
      eC = findall(x->elabels[x]==i, 1:nedges2)
      edgecnt[i] = length(findall(x->elabels[x]==i, 1:nedges2))
      push!(clusters, C)
      edgesi = random_hypergraph(C, edgecnt[i], hep2; kwargs...)
      append!(edges, edgesi)
    end
    return unique!(sort!.(edges))
  end
  
  #n = 10000
  #Hlist = hypergraph_SBM(n, 100, 10000, 50000, 0.8, 0.8)

#=
Test about Jamie Haddock's idea for opinion dynamics with 
  some type of nonlin
=#  


function project(edges)
  gedges = Tuple{Int,Int}[] 
  for e in edges
    for u in e
      for v in e
        if u != v
          push!(gedges, (u,v))
        end
      end
    end
  end
  return gedges 
end 

function edges2sparse(edges)
  n = maximum(e->max(e[1], e[2]), edges)
  A = sparse(first.(edges), last.(edges), 1, n, n)
  return A 
end 

## code for pairwise spatial model
using NearestNeighbors, Distributions, SparseArrays
function spatial_graph_edges(n::Integer,d::Integer;degreedist=LogNormal(log(4),1))
  xy = rand(d,n)
  T = BallTree(xy)
  # form the edges for sparse
  ei = Int[]
  ej = Int[]
  for i=1:n
    deg = min(ceil(Int,rand(degreedist)),n-1)
    idxs = knn(T, xy[:,i], deg+1)[1]
    for j in idxs
      if i != j
        push!(ei,i)
        push!(ej,j)
      end
    end
  end
  return xy, ei, ej
end
function spatial_network(n::Integer, d::Integer; degreedist=LogNormal(log(3),1))
  xy, ei, ej = spatial_graph_edges(n, d;degreedist=degreedist)
  A = sparse(ei,ej,1,n,n)
  return max.(A,A'), xy
end

using CairoMakie
function plotgraph(A::SparseMatrixCSC,xy::AbstractArray{T,2};kwargs...) where T
  px,py = zeros(T,0),zeros(T,0)
  P = [px,py]
  rows = rowvals(A)
  skip = NaN.*xy[:,begin] # first row
  for j=1:size(A,2) # for each column
    for nzi in nzrange(A, j)
      i = rows[nzi]
      if i > j
        push!.(P, @view xy[:,i])
        push!.(P, @view xy[:,j])
        push!.(P, skip)
      end
    end
  end
  plot(px,py;framestyle=:none,legend=false,kwargs...)
end

## hypergraph spatial model 
function spatial_hypergraph_edges(n::Integer,d::Integer;degreedist=LogNormal(log(4),1))
  xy = rand(d,n)
  T = BallTree(xy)
  # form the edges for sparse
  edges = Vector{Int}[]
  for i=1:n
    deg = min(ceil(Int,rand(degreedist)),n-1)
    # now we need to sample sizes for each edge... 

    idxs = knn(T, xy[:,i], deg+1)[1]

    # want a random subset of these indices...
    for j in idxs
      if i != j
        e = [ei,ej]
        while rand() < pedge 
          j2 = i # start in the 'repeat' loop... 
          while j2 == i 
            j2 = rand(indx)
            push!(e, j2)
          end 
        end
      end
    end
  end
  return xy, ei, ej
end


## File I/O
function save_hedges(hedges,fpath)
  if !endswith(fpath,".txt")
    fpath = "$fpath.txt"
  end
  open(fpath, "w") do io
    for edge in hedges
      writedlm(io, [edge], ',')
    end
  end
end
function read_hedges(filename::String)
  result = Vector{Vector{Int}}()
  open(filename, "r") do io
      for line in eachline(io)
          # Split the line by commas and parse each element to Int
          inner_vec = parse.(Int, split(line, ","))
          push!(result, inner_vec)
      end
  end
  return result
end
