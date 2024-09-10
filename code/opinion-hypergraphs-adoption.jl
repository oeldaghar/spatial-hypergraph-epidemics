##
# generate a random hypergraph
using Random, MatrixNetworks
include("hypergraph-tools.jl")
include("spatial-hypergraph.jl")
using Random

Random.seed!(0)
n = 50
hedges, xy = spatial_hypergraph_edges(n, 2)
scatter(xy[1,:],xy[2,:])

segs, meanpts = graphlines_spatial_hypergraph(hedges, xy; degreedist = 8:8)
f = linesegments(segs)
scatter!(xy[1,:], xy[2,:])
scatter!(meanpts, markersize = 15)
f
#heatmap(G)
##
function _projection!(x, e) 
  # get the average of xi  
  n = length(e) 
  xmean = sum(i->x[i]/n, e)
  #xmean = xmean/abs(xmean) 
  for i in e
    xi = x[i] 
    x[i] = xi*xmean
    #x[i] = x[i]/abs(x[i])
    #x[i] = (xmean + xi)/2
  end 
  return x
end

function opinion_step!(x, H)
  j = rand(1:length(H)) # find a random edge
  e = H[j] 
  # compute the projection step 
  _projection!(x, e) 
end 


Random.seed!(4)
#x = exp(pi/2*im)*ones(n)
x = 0*ones(n)
x[(xy[1,:].<0.2) .& (xy[2,:] .<0.2)] .= 1
function run_steps(x, H, nsteps=100)
  X = zeros(eltype(x), length(x), nsteps)
  @show typeof(X)
  for i in 1:nsteps
    copyto!(@view(X[:,i]), opinion_step!(x, H))
  end 
  return X 
end 
X = run_steps(x, hedges)
T = X



##
fig = Figure(size=(1024,1024))
ax = Axis(fig[1, 1])

sl_i = Slider(fig[2, 1], range = 1:size(T,2), startvalue = 1)

colors = lift(sl_i.value) do i
    T[:,i]
end

linesegments!(segs, linewidth=0.5)
scatter!(xy[1,:], xy[2,:], color=colors)
fig
