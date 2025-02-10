"""
script for showing differences in the embedding X

uniform embedding
mixture of guassianas 
"""

using Distributed
# Distributed.addprocs(10)
@everywhere include("spatial-hypergraph.jl")
@everywhere using Distributions, Random
@everywhere using StatsBase
@everywhere using ProgressMeter
@everywhere using MatrixNetworks, Graphs, SparseArrays
@everywhere using SparseArrays

# Plotting 
using Plots, LaTeXStrings, Colors
using Measures


# what do we want to measure?
# hypergraph size distribution? - edge counts 


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

function quantile_figure(data,alphas)
    quantiles = [0.1,0.25,0.5,0.75,0.9]
    linewidths = [0.5, 1, 2, 1, 0.5]
    colors = [:lightgrey, :grey, :black, :grey, :lightgrey]

    # total edges 
    f = Plots.plot(xlabel=L"\alpha", ylabel="Total hyperedges", 
                framestyle=:box,
                thickness_scaling=1.2,
                guidefontsize=14,
                tickfontsize=12,
                tickdirection=:out)
    for (q, lw, c) in zip(quantiles, linewidths, colors)
        nedges_q = quantile.(eachrow(data), q)
        Plots.plot!(f, alphas, nedges_q, label="", linewidth=lw, color=c,    
            yscale=:log10)
    end
    return f
end


function _remove_tick_labels(f)
    curr_yticks = Plots.yticks(f)
    new_yticks = (curr_yticks[1],["" for i=1:lastindex(curr_yticks[1])])
    Plots.plot!(f,yticks=new_yticks,ylabel="")
end    


alphas = range(0,2,25)
n = 10000
d = 10
ntrials = 25
figs = []
for alpha_func in [get_func,alpha_func_linear,alpha_func_no_dim]
    trials = _get_plotting_data(n,d,alphas,ntrials,alpha_func)
    nedges = reduce(hcat,trials)
    push!(figs,quantile_figure(nedges,alphas))
end

# put figure together 
plt = Plots.plot(figs...,layout = (1,3))
# touch up margins and label 
Plots.plot!(plt,plot_title = "n=$n, d=$d",link=:all,size=(1800,500),
                plot_titlefontsize=22)
_remove_tick_labels(plt[2])
_remove_tick_labels(plt[3])
Plots.plot!(plt,bottom_margin=8mm)
Plots.plot!(plt[1],left_margin=8mm)
Plots.plot!(plt,top_margin=8mm,dpi = 1000)

Plots.savefig(plt,"data/output/figures/final/alpha-figure.pdf")
Plots.savefig(plt,"data/output/figures/final/alpha-figure.png")