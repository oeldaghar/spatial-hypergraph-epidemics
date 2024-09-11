#load code on all workers
using Distributed
using Statistics
addprocs(2)

println("LOADING CODE FILES AND PACKAGES")
@everywhere include("sirs.jl")

fnames = filter(x->endswith(x,".txt"),readdir("data/hypergraphs/"))

fname = fnames[1]
hedges = read_hedges(joinpath("data/hypergraphs/$(fnames[1])"))
edges = project(hedges)


nnodes = maximum(x->maximum(x),edges)
seed_nodes = Distributions.sample(1:nnodes,100,replace=false)

pairwise_results = parallel_sirs(edges,seed_nodes)
hyper_results = parallel_sirs_hyper(hedges,seed_nodes)

pairwise_data = hcat(pairwise_results...)';
hyper_data = hcat(hyper_results...)';

function make_simulation_bands(pairwise_data,hyper_data)
    #figure layout 
    #put figures side by side 
    fig = Figure(resolution=(1000, 400))
    left_column = fig[1, 1] = GridLayout()
    right_column = fig[1, 2] = GridLayout()
    # Create axes in each column
    ax1 = Axis(left_column[1, 1], 
                title="HigherOrder Model",
                xlabel="Time",
                ylabel="Total Infections")
    ax2 = Axis(right_column[1, 1],
                title="Pairwise Model",
                xlabel="Time",
                ylabel="Total Infections")
    # Make plotting data 
    #get average and 2.5-97.5 quantiles 
    pairwise_avgs = mapslices(x->mean(x),pairwise_data,dims=1)
    pairwise_ylow = mapslices(x->quantile(x,0.01),pairwise_data,dims=1)
    pairwise_yhigh = mapslices(x->quantile(x,0.99),pairwise_data,dims=1)
    pairwise_avgs,pairwise_ylow,pairwise_yhigh = vec(pairwise_avgs),vec(pairwise_ylow),vec(pairwise_yhigh)

    #hyper data 
    hyper_avgs = mapslices(x->mean(x),hyper_data,dims=1)
    hyper_ylow = mapslices(x->quantile(x,0.01),hyper_data,dims=1)
    hyper_yhigh = mapslices(x->quantile(x,0.99),hyper_data,dims=1)
    hyper_avgs,hyper_ylow,hyper_yhigh = vec(hyper_avgs),vec(hyper_ylow),vec(hyper_yhigh)
    
    #render 
    lines!(ax1, hyper_avgs,color=(:blue,1.0))
    band!(ax1, 1:lastindex(hyper_avgs), hyper_ylow, hyper_yhigh, color=(:blue, 0.25))
    
    lines!(ax2, pairwise_avgs,color=(:red,1.0))
    band!(ax2, 1:lastindex(pairwise_avgs), pairwise_ylow, pairwise_yhigh, color=(:red, 0.25))

    linkaxes!(ax1, ax2)
    # Display the figure
    return fig,ax1,ax2     
end


fig,ax1,ax2 = make_simulation_bands(pairwise_data,hyper_data)
fig

xlims!(ax1,100,5000)
xlims!(ax2,100,5000)
ylims!(ax1,0,500)
ylims!(ax2,0,500)
display(fig)
