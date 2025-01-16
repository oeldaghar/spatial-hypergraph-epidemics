using Distributed
using JSON 
# addprocs(2)

# include stuff
@everywhere include("../data-io.jl")
@everywhere include("../hypergraph-tools.jl")
@everywhere include("../spatial-hypergraph.jl")
@everywhere include("../sirs.jl")

# test with smaller graph 
gname = "spatial-hypergraph-5000-2"
fnames = filter!(x->endswith(x,".txt"), get_fnames(gname))
nnodes = parse(Int,split(fnames[1],'-')[3])

# need to specify the following parameters 
# seed_nodes: nodes from which to seed the epidemic
# beta: pairwise infection probability to use  
# gamma: node recovery probability
# delta: node relapse probability (R->S) 
# exo: term for exogenous infections 
# tmax: total number of timesteps 
# hyper_beta_func: normalization term to use for hyperedges 

## ALL PARAMETERS 
nsamples_per_seeds = 10
initial_infected_fractions = [0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 0.95]
# epidemic parameters 
beta = vcat(1e-3:1e-3:1e-2, 2e-2:1e-2:1e-1, 2e-1:1e-1:9e-1)
gamma,delta,exo = [5e-2],[1/20],[5/1e6]
tmax = [365*10]
hyper_beta_func = ["linear","sqrt"]

# PARAMETERS FOR TESTING
# nsamples_per_seeds = 1
# initial_infected_fractions = [0.05, 0.1]
# # epidemic parameters 
# beta = [1e-2] #vcat(1e-3:1e-3:1e-2, 2e-2:1e-2:1e-1, 2e-1:1e-1:9e-1)
# gamma,delta,exo = [5e-2],[5e-2],[0.0]
# tmax = [365*10]
# hyper_beta_func = ["linear"]

# random seed for choosing seed nodes for infection
rseed=abs(rand(Int))
Random.seed!(rseed)
# seed nodes. same seeds as we vary alpha.
num_nodes_sampled = round.(Int,nnodes*initial_infected_fractions)
initial_seeds = Vector{Vector{Int}}()
for nseeds in num_nodes_sampled
    for i=1:nsamples_per_seeds    
        seed_nodes = Distributions.sample(1:nnodes,nseeds,replace=false)
        push!(initial_seeds,seed_nodes)
    end
end
# compartmental model parameters 

parallel_time = 0

# MAIN LOOP
for (i,fname) in enumerate(fnames)
    println("WORKING ON GRAPH $i OF $(length(fnames))")
    hedges = read_hedges(fname)

    println("RUNNING PARALLEL EPIDEMIC CODE...")
    start_time = time()
    # main parallel function call
    hyper_results = parallel_sirs_hyper(hedges,initial_seeds;
                                            beta=beta, gamma=gamma, delta=delta, exo=exo, 
                                            tmax=tmax, hyper_beta_func=hyper_beta_func)
    end_time = time()
    delta_time = end_time-start_time
    global parallel_time += delta_time
    println("GRAPH $i ELAPSED TIME USING $(nworkers()) WORKERS: $(delta_time) seconds")
    infection_information = vcat(hyper_results[1]...)
    epi_params = vcat(hyper_results[2]...)

    data_path = "data/hysteresis/sirs/"
    dst_file = "$(splitext(basename(fname))[1]).json"
    dst_path = joinpath(data_path,dst_file)
    # keep files for now, later move to just overwrite with a flag 
    counter = 1
    while isfile(dst_path)
        dst_file = "$(splitext(basename(fname))[1])-counter_$counter.json"
        dst_path = joinpath(data_path,dst_file)
        counter+=1
    end

    if !(ispath(data_path))
        mkpath(data_path)
    end

    # write to JSON 
    open(dst_path, "w") do io
        output = Dict(
            # epidemic parameters 
            "initial_num_infected_nodes" => map(x->lastindex(x[1]),epi_params), #first.(epi_params),
            "beta" => map(x->x[2],epi_params),
            "gamma" => map(x->x[3],epi_params),
            "delta" => map(x->x[4],epi_params),
            "exo" => map(x->x[5],epi_params),
            "tmax" => map(x->x[6],epi_params),
            "hyper_beta_func" => map(x->x[7],epi_params),
            # infection information 
            "infections" => first.(infection_information),
            "ninfected_hyperedge" => map(x->x[2],infection_information),
            "ntransmissions_hyperedge" => map(x->x[3],infection_information),
        )
        JSON.print(io, output)
    end
end 
println("TOTAL PARALLEL TIME: $parallel_time.\nTOTAL WORKERS: $(nworkers())")

# TBW
# want one diagram to be alpha vs beta vs infected fraction 
# need to find transition or use a granular enough parameter regime

# try to read it back in before running lots of experiments 
# length.(JSON.parsefile(dst_path)["seed_nodes"])
# fpath = "data/hysteresis/sirs/spatial-hypergraph-5000-2-0.0-1-counter_1-TEST.json"
# # fpath = "data/hysteresis/sirs/spatial-hypergraph-5000-2-0.0-1.json"
# data = JSON.parsefile(fpath)
# data["infections"]

# data["infections"]
# hcat(data["ninfected_hyperedge"]...)
