## Need a spatial graph model
using Clustering 
using NearestNeighbors

# This function is designed to interpolate between pure hypergraph at 2
# and pure graph at 0.
#get_func(alpha) = (d,deg) -> (d/2 - d/sqrt(deg))*alpha^2 + (2*d/sqrt(deg) - d/2)*alpha
#  replace sqrt(deg) w/ deg^(1/d)
function get_func(alpha) 
  if alpha <= 1
    return (dist,deg,dim) -> (alpha^(1/dim)*dist)/(deg)^(1/(2*dim))
  else
    return (dist,deg,dim) -> dist/(deg)^(1/(2*dim)) + (alpha-1)*(dist-dist/(deg)^(1/(2*dim)))
  end
end 

function hypergraph_edges(X;degreedist=LogNormal(log(3),1),radfunc=(dist,deg,dim) -> dist/(deg)^(1/(2*dim)))
  T = BallTree(X)
  # form the edges for sparse
  edges = Vector{Int}[]
  n = lastindex(X,2)
  l = lastindex(X,1)
  for i=1:n #X
    deg = min(ceil(Int,rand(degreedist)),n-1)
    idxs, dists = knn(T, X[:,i], deg+1)
    if deg > 1 
      maxdist = maximum(dists) 
      pts = @view X[:,idxs]
      rad = radfunc(maxdist,deg,l)
      # if !(rad ≈ maxdist)
      #   println("rad: ", rad, " maxdist: ", maxdist)
      # end
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

function hypergraph_edges(X;degs::Vector{S}=rand(LogNormal(log(3),1),lastindex(X,2)),radfunc=(dist,deg,dim) -> dist/(deg)^(1/(2*dim))) where S
  T = BallTree(X)
  # form the edges for sparse
  edges = Vector{Int}[]
  n = lastindex(X,2)
  l = lastindex(X,1)
  for i=1:n #X
    deg = degs[i]
    idxs, dists = knn(T, X[:,i], deg+1)
    if deg > 1 
      maxdist = maximum(dists) 
      pts = @view X[:,idxs]
      rad = radfunc(maxdist,deg,l)
      # if !(rad ≈ maxdist)
      #   println("rad: ", rad, " maxdist: ", maxdist)
      # end
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

# hypergraph_edges(X, deg_list::Vector{S}; radfunc) where S = hypergraph_edges(X, degs=deg_list, radfunc=radfunc)

function spatial_hypergraph_edges(n::Integer,d::Integer;kwargs...)
  X = rand(d,n)
  return hypergraph_edges(X;kwargs...)
end

function hypergraph_edges(X,deg_list::Vector{S};radfunc=(dist,deg,dim) -> dist/(deg)^(1/(2*dim))) where S
  T = BallTree(X)
  edges = Vector{Int}[]
  dim = size(X,1)
  for i in axes(X,2)
    # deg = min(ceil(Int,rand(degreedist)),n-1)
    deg = deg_list[i]
    idxs, dists = knn(T, X[:,i], deg+1)
    if deg > 1 
      maxdist = maximum(dists) 
      pts = @view X[:,idxs]
      rad = radfunc(maxdist,deg,dim)
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
    end 
  end 
  return edges
end

#example usage 
# alphas = [0,2.0]
# X = rand(2,5000) #d = dimension, n = num nodes  
# hypergraph_edges(X;radfunc=get_func(alphas[1]))

##
# using CairoMakie
# function graphlines_spatial_hypergraph(hedges, xy::AbstractMatrix{T};kwargs...) where T
#   #px,py = zeros(T,0),zeros(T,0)
#   #P = [px,py]
#   #skip = NaN.*xy[:,begin] # first col
#   segs = Tuple{Point2f,Point2f}[] 
#   meanpts = Point2f[] 
#   for e in hedges
#     if length(e) == 2
#       #push!.(P, @view xy[:,e[1]])
#       #push!.(P, @view xy[:,e[2]])
#       #push!.(P, skip)
#       push!(segs, (Point2f(xy[:,e[1]]), Point2f(xy[:,e[2]])))
#     else
#       # find the centroid... 
#       hsize = length(e) 
#       mx = sum( i->xy[1, i]/hsize , e)
#       my = sum( i->xy[2, i]/hsize , e)
#       mp = Point2f(mx,my) 
#       push!(meanpts, mp)
#       for v in e
#         push!(segs, (mp, Point2f(xy[:,v])))
#       end 
#     end
#   end 
#   return segs, meanpts
# end
