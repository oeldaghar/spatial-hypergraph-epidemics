using Distributed
using JSON 
using JSON3

# include auxillary files 
@everywhere include("../data-io.jl")
@everywhere include("../hypergraph-tools.jl")
@everywhere include("../spatial-hypergraph.jl")
@everywhere include("../sirs.jl")

# where files are stored 
DATA_PATH = "data/hysteresis/sirs/scratch-v3/"
if !(ispath(DATA_PATH))
    mkpath(DATA_PATH)
end

# need to specify the following parameters 
# seed_nodes: nodes from which to seed the epidemic
# beta: pairwise infection probability to use  
# gamma: node recovery probability
# delta: node relapse probability (R->S) 
# exo: term for exogenous infections 
# tmax: total number of timesteps 
# hyper_beta_func: normalization term to use for hyperedges 

## ALL PARAMETERS 
nsamples_per_seeds = 5
# initial_infected_fractions = [0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 0.95]
initial_infected_fractions = [0.05, 0.1]#, 0.15] #, 0.2, 0.25]
# epidemic parameters 
# beta = vcat(1e-3:1e-3:1e-2, 2e-2:1e-2:1e-1, 2e-1:1e-1:9e-1)
beta = [1e-2, 2e-2, 3e-2, 4e-2, 5e-2 ,0.1, 0.9]
# gamma,delta,exo = [5e-2],[1/20],[5/1e6]
gamma, delta, exo = [0.05], [0.01, 0.02, 0.03, 0.04], [5/1e6] # from the conference paper
# tmax = [365*10]
tmax = [2000]
# hyper_beta_func = ["linear","sqrt","pairwise"] # full set but running incremental diffusions
hyper_beta_func = ["pairwise", "linear", "sqrt"]

# PARAMETERS FOR TESTING
# nsamples_per_seeds = 1
# initial_infected_fractions = [0.05, 0.1]
# # epidemic parameters 
# beta = [1e-2,3e-2,5e-2] #vcat(1e-3:1e-3:1e-2, 2e-2:1e-2:1e-1, 2e-1:1e-1:9e-1)
# gamma,delta,exo = [5e-2],[5e-2],[0.0]
# tmax = [365*10]
# hyper_beta_func = ["linear","sqrt","pairwise"]

# Function to write a single epidemics result to the JSONL file
function write_result(io, params, info)
    result = Dict(
        "initial_num_infected_nodes" => lastindex(params[1]),
        "beta" => params[2],
        "gamma" => params[3],
        "delta" => params[4],
        "exo" => params[5],
        "tmax" => params[6],
        "hyper_beta_func" => params[7],
        "infections" => info[1],
        "ninfected_hyperedge" => info[2],
        "ntransmissions_hyperedge" => info[3]
    )
    JSON3.write(io, result)
    write(io, '\n')  # Add a newline after each JSON object
end

gnames = [
        #smaller graphs         
        "spatial-hypergraph-5000-2", 
        "spatial-hypergraph-5000-5",
        #larger graphs 
        "spatial-hypergraph-50000-2", 
        "spatial-hypergraph-50000-5",
]

summary_info = []
# GRAPH LOOP
for gname in gnames 
    fnames = filter!(x->endswith(x,".txt"), get_fnames(gname))
    fnames = filter(x->occursin("newalpha",x),fnames)
    nnodes = parse(Int,split(fnames[1],'-')[3])

    # random seed for choosing initial infected nodes
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

    parallel_time = 0
    # MAIN LOOP
    for (i,fname) in enumerate(fnames)
        other_start_time = time()
        println("WORKING ON GRAPH $i OF $(length(fnames))")
        hedges = read_hedges(fname)

        println("RUNNING PARALLEL EPIDEMIC CODE...")
        parallel_start_time = time()
        # MAIN PARALLEL FUNCTION CALL 
        hyper_results = parallel_sirs_hyper(hedges,initial_seeds;
                                                beta=beta, gamma=gamma, delta=delta, exo=exo, 
                                                tmax=tmax, hyper_beta_func=hyper_beta_func)
        parallel_end_time = time()
        delta_parallel_time = parallel_end_time - parallel_start_time
        parallel_time += delta_parallel_time
        println("GRAPH $i ELAPSED TIME USING $(nworkers()) WORKERS: $(delta_parallel_time) seconds")

        # TODO test writing this to another format like Parquet 
        dst_fname = "$(splitext(basename(fname))[1]).jsonl"    
        jsonl_file = joinpath(DATA_PATH, dst_fname)
        println("SAVING RESULTS TO FILE: $jsonl_file")
        stime = time()
        open(jsonl_file, "a") do io
            # infection_information = reduce(vcat,hyper_results[1])
            # epi_params = reduce(vcat,hyper_results[2])
            for (params, info) in zip(reduce(vcat,hyper_results[2]), reduce(vcat,hyper_results[1]))
                write_result(io, params, info)
            end
        end
        etime = time()
        println("WROTE TO FILE IN $(etime-stime) SECONDS")
    end 
    println("TOTAL PARALLEL TIME: $parallel_time.\nTOTAL WORKERS: $(nworkers())")
    push!(summary_info,parallel_time)
end

for (ind,ptime) in enumerate(summary_info)
    println("TOTAL PARALLEL TIME FOR $(gnames[ind]): $ptime SECONDS")
end
