using Distributed
addprocs(17)

@everywhere include("hypergraph-tools.jl")
@everywhere include("spatial-hypergraph.jl")

d = 2
n = 10000
using Random
Random.seed!(1234)
X = rand(d,n) #d = dimension, n = num nodes  
degreedist=LogNormal(log(3),1)
degs = rand(degreedist,n)
degs = min.(ceil.(Int,degs),n-1)
# This function is designed to interpolate between pure hypergraph at 2
# and pure graph at 0.
#get_func(alpha) = (d,deg) -> (d/2 - d/sqrt(deg))*alpha^2 + (2*d/sqrt(deg) - d/2)*alpha
@everywhere function get_func(alpha) 
  if alpha <= 1
    return (d,deg) -> alpha*d/sqrt(deg)
  else
    return (d,deg) -> d/sqrt(deg) + (alpha-1)*(d-d/sqrt(deg))
  end
end 
alphas = range(0,2,length=25)
#alphas = [2.0]
@everywhere X,alphas = $X,$alphas
graphs = pmap(alpha -> hypergraph_edges(X;degs=degs,radfunc=get_func(alpha))[1], alphas)

##
nedges = map(x->length(x), graphs)

##
using Plots
Plots.plot(alphas, nedges, label="Number of edges", xlabel="alpha", ylabel="Number of edges", legend=:topleft)

## Now let's get this for multiple realizations...
n = 50000
d = 2 
alphas = range(0,2,length=25)
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

## Now let's plot the average number of edges...
function plot_data(trials)
  nedges = reduce(hcat, trials) # now we have alpha as row index and trial as column index
  # I want to take the quantiles over the columns... 
  quantiles = [0.1,0.25,0.5,0.75,0.9]
  linewidths = [0.5, 1, 2, 1, 0.5]
  colors = [:lightgrey, :grey, :black, :grey, :lightgrey]

  p = Plots.plot(xlabel="α", ylabel="Number of edges")
  for (q, lw, c) in zip(quantiles, linewidths, colors)
    nedges_q = quantile.(eachrow(nedges), q)
    Plots.plot!(p, alphas, nedges_q, label="", linewidth=lw, color=c)
  end
  p
end
plot_data(trials)

## Now I want to project each graph
n = 50000
d = 2 
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

##
projected_trials_squared = @showprogress map(x->run_projected_trial(scalefunc=x->1/x^2), 1:25)  
p1 = plot_data(projected_trials_squared)
Plots.plot!(p1, ylabel="Weighted projected edge volume")
p1


##
projected_trials_sqrt = @showprogress map(x->run_projected_trial(scalefunc=x->1/sqrt(x)), 1:25)  
p2 = plot_data(projected_trials_sqrt)
Plots.plot!(p2, ylabel="Weighted projected edge volume")
p2


##
projected_trials = @showprogress map(x->run_projected_trial(), 1:25)  
p3 = plot_data(projected_trials)
Plots.plot!(p3, ylabel="Weighted projected edge volume")
p3

##


using Measures
#making combined plot
figs = Plots.plot(p3,p2,p1,
            layout=(1,3),size=(1400,400),
            bottom_margin=8Measures.mm,left_margin=8Measures.mm,
)
Plots.savefig(figs,"data/output/figures/final/projected-degree.pdf")
Plots.savefig(figs,"data/output/figures/final/projected-degree.png")



# ##figures for degree distribution by size 
# function run_projected_trial(;scalefunc=x->1/x)
#     X = rand(d,n)
#     degreedist=LogNormal(log(3),1)
#     degs = rand(degreedist,n)
#     degs = min.(ceil.(Int,degs),n-1)

#     graphs = pmap(alpha -> hypergraph_edges(X;degs=degs,radfunc=get_func(alpha))[1], alphas)
#     # map(x->,graphs)
#     nedges = map(x->_project_graph_avgdegrees(x;scalefunc), graphs)
#     return nedges   
# end 


manual_xlims = (-1,130)
p1 = Plots.histogram(length.(graphs[1]),
            normalize=:probability,label=false,
            xlabel="Hyperedge Size",
            ylabel="Fraction of Hyperedges")
Plots.xlims!(manual_xlims)
p2 = Plots.histogram(length.(graphs[13]),
            normalize=:probability,label=false,
            xlabel="Hyperedge Size",
            ylabel="Fraction of Hyperedges")
Plots.xlims!(manual_xlims)
p3 = Plots.histogram(length.(graphs[end]),
            normalize=:probability,label=false,
            xlabel="Hyperedge Size",
            ylabel="Fraction of Hyperedges")
Plots.xlims!(manual_xlims)

Plots.plot(p1,p2,p3,link=:both,layout=(1,3),size=(1400,400),
            bottom_margin=8Measures.mm,
            left_margin=8Measures.mm)

function logscale_histogram(xs,ymin=1/length(xs),normalize_y::Bool=true;logx::Bool=false,maxbins::Int=5000)
    maxbins = min(maxbins,ceil(1+max(abs.(extrema(xs))...)))
    h = float(StatsBase.fit(Histogram,xs,nbins=maxbins))
    if normalize_y
        h.weights./=length(xs)
        ymax = 1
    else
        ymax = length(xs)
    end

    #a little padding to make plots look better
    bins = [collect(h.edges[1]);maximum(h.edges[1])+step(h.edges[1])]
    w = [h.weights;0] .+ 0.5*ymin
    if logx
        f = Plots.plot(bins,w,seriestype=:barbins,label="",
                xlims=(1,2*maximum(xs)),xscale=:log10,
                ylims=(ymin,ymax),yscale=:log10,bar_width=step(h.edges[1]),lims=:round)
    else
        f = Plots.plot(bins,w,seriestype=:barbins,label="",
                xlims=(0,maximum(xs)+2),ylims=(ymin,ymax),yscale=:log10,lims=:round)
    end
    return f
end
using StatsBase
p1 = logscale_histogram(length.(graphs[1]))
p2 = logscale_histogram(length.(graphs[13]))
p3 = logscale_histogram(length.(graphs[end]))


counts(length.(graphs[end]))