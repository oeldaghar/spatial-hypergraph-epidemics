#load code on all workers
using Distributed
using Statistics
using DelimitedFiles
addprocs(5)

println("LOADING CODE FILES AND PACKAGES")
@everywhere include("sirs.jl")
include("data-io.jl")

fnames = filter!(x->endswith(x,".txt"),readdir("data/hypergraphs/"))
fnames = get_fnames(get_gname(fnames[1]))

#sample random nodes across all graphs
n = parse(Int,split(fnames[1],'-')[3])
num_nodes_sampled = 10
Random.seed!(42)
seed_nodes = Distributions.sample(1:n,num_nodes_sampled,replace=false)

beta,gamma,delta,exo,tmax = 0.9,1e-1,1/100,5/1e6,365*10
# hyperr_beta_func = "linear" # "sqrt", "squared"
# run parallel code and save data 
for hyper_beta_func in ["linear","sqrt","squared"]
    @showprogress for fname in fnames 
        #load data 
        hedges = read_hedges(fname)
        println("RUNNING PARALLEL EPIDEMIC CODE...")
        hyper_results = parallel_sirs_hyper(hedges,seed_nodes;
                                                beta=beta, gamma=gamma, delta=delta, exo=exo, 
                                                tmax=tmax, hyper_beta_func=hyper_beta_func)
        hyper_data = hcat(hyper_results...)'
        
        gname = split(fname,'/')[end]
        if endswith(gname,".txt")
            gname = gname[1:end-4]
        end 
        file_name = "data/output/sirs/$gname%sirs-$beta-$gamma-$delta-$exo-$tmax-$hyper_beta_func.txt"
        open(file_name, "w") do io
            writedlm(io, hyper_data)
        end
    end
end



# make plot

# hyper_results = parallel_sirs_hyper(hedges,seed_nodes)
# hyper_data = hcat(hyper_results...)';

# function make_simulation_bands(hyper_data)
#     #figure layout 
#     #put figures side by side 
#     fig = Figure()
#     # Create axes in each column
#     ax = Axis(fig[1,1], 
#                 title="HigherOrder Model",
#                 xlabel="Time",
#                 ylabel="Total Infections")
#     # Make plotting data 
#     hyper_avgs = mapslices(x->mean(x),hyper_data,dims=1)
#     hyper_ylow = mapslices(x->quantile(x,0.01),hyper_data,dims=1)
#     hyper_yhigh = mapslices(x->quantile(x,0.99),hyper_data,dims=1)
#     hyper_avgs,hyper_ylow,hyper_yhigh = vec(hyper_avgs),vec(hyper_ylow),vec(hyper_yhigh)
    
#     #render 
#     lines!(ax, hyper_avgs,color=(:blue,1.0))
#     band!(ax, 1:lastindex(hyper_avgs), hyper_ylow, hyper_yhigh, color=(:blue, 0.25))
    
#     # Display the figure
#     return fig,ax
# end


# fig,ax = make_simulation_bands(hyper_data)
# fig

# xlims!(ax1,100,5000)
# xlims!(ax2,100,5000)
# ylims!(ax1,0,500)
# ylims!(ax2,0,500)
# display(fig)
