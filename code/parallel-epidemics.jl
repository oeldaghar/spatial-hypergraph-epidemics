#load code on all workers
using Distributed
using Statistics
using DelimitedFiles
addprocs(25)

println("LOADING CODE FILES AND PACKAGES")
@everywhere include("sirs.jl")
include("data-io.jl")

#parameters to tune here 
num_nodes_sampled = 25
beta,gamma,delta,exo,tmax = 0.9,1e-1,1/100,5/1e6,365*10

for gname in ["spatial-hypergraph-5000-2","spatial-hypergraph-50000-2"]
    fnames = get_fnames(gname)
    #specify graph name to use 
    fnames = filter!(x->endswith(x,".txt"),fnames)
    # run parallel code and save data 
    for hyper_beta_func in ["linear","sqrt","squared"]
        @showprogress for fname in fnames 
            #load data 
            rseed=abs(rand(Int))
            Random.seed!(rseed)
            n = parse(Int,split(fnames[1],'-')[3])
            seed_nodes = Distributions.sample(1:n,num_nodes_sampled,replace=false)

            hedges = read_hedges(fname)
            println("RUNNING PARALLEL EPIDEMIC CODE...")
            hyper_results = parallel_sirs_hyper(hedges,seed_nodes;
                                                    beta=beta, gamma=gamma, delta=delta, exo=exo, 
                                                    tmax=tmax, hyper_beta_func=hyper_beta_func)
            hyper_data = hcat(hyper_results...)'
            println("SAVING RESULTS...")
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
end