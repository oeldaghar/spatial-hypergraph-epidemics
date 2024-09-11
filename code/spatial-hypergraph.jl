## Need a spatial graph model
using Clustering 
function spatial_hypergraph_edges(n::Integer,d::Integer;degreedist=LogNormal(log(3),1))
  xy = rand(d,n)
  T = BallTree(xy)
  # form the edges for sparse
  edges = Vector{Int}[]
  for i=1:n
    deg = min(ceil(Int,rand(degreedist)),n-1)
    idxs, dists = knn(T, xy[:,i], deg+1)
    if deg > 1 
      maxdist = maximum(dists) 
      pts = @view xy[:,idxs]
      rad = maxdist/sqrt(deg)
      clusters = dbscan(pts, rad)
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
  return edges, xy
end


##
using CairoMakie
function graphlines_spatial_hypergraph(hedges, xy::AbstractMatrix{T};kwargs...) where T
  #px,py = zeros(T,0),zeros(T,0)
  #P = [px,py]
  #skip = NaN.*xy[:,begin] # first col
  segs = Tuple{Point2f,Point2f}[] 
  meanpts = Point2f[] 
  for e in hedges
    if length(e) == 2
      #push!.(P, @view xy[:,e[1]])
      #push!.(P, @view xy[:,e[2]])
      #push!.(P, skip)
      push!(segs, (Point2f(xy[:,e[1]]), Point2f(xy[:,e[2]])))
    else
      # find the centroid... 
      hsize = length(e) 
      mx = sum( i->xy[1, i]/hsize , e)
      my = sum( i->xy[2, i]/hsize , e)
      mp = Point2f(mx,my) 
      push!(meanpts, mp)
      for v in e
        push!(segs, (mp, Point2f(xy[:,v])))
      end 
    end
  end 
  return segs, meanpts
end
