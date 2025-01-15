include("data-io.jl")
include("hypergraph-tools.jl")
include("spatial-hypergraph.jl")
include("sirs.jl")

gname = "spatial-hypergraph-5000-2"

fnames = get_fnames(gname)
#specify graph name to use 
fnames = filter!(x->endswith(x,".txt"),fnames)

#load data 
rseed=abs(rand(Int))
Random.seed!(rseed)
n = parse(Int,split(fnames[1],'-')[3])

#parameters to test 
# seed_nodes::Vector{Vector{S}} = [[1]],
# beta::Vector{T} = [5e-2], 
# gamma::Vector{T} = [5e-2], 
# delta::Vector{T} = [1/365], 
# exo::Vector{T} = [5/100000], 
# tmax::Vector{S} = [365*10], 
# hyper_beta_func::Vector{String} = ["linear"]) where {T, S}

nsamples_per_seeds = 10
initial_infected_fractions = [0.1,0.75]
num_nodes_sampled = round.(Int,n*initial_infected_fractions)
initial_seeds = Vector{Vector{Int}}()
for nseeds in num_nodes_sampled
    for i=1:nsamples_per_seeds    
        seed_nodes = Distributions.sample(1:n,nseeds,replace=false)
        push!(initial_seeds,seed_nodes)
    end
end
beta,gamma,delta,exo,tmax = [0.01],[5e-2],[1/20],[5/1e6],[365*10]
hyper_beta_func = ["linear","sqrt"]

fname = fnames[end]
hedges = read_hedges(fname)

println("RUNNING PARALLEL EPIDEMIC CODE...")
hyper_results = parallel_sirs_hyper(hedges,initial_seeds;
                                        beta=beta, gamma=gamma, delta=delta, exo=exo, 
                                        tmax=tmax, hyper_beta_func=hyper_beta_func)

hyper_data = hcat(hyper_results[1]...)'
epi_params = hyper_results[2]

# average of last 1000 time steps 
rolling_avg = mapslices(mean,hyper_data[:,end-1000:end],dims=2)

a,b = extrema(rolling_avg)
println("minimum : $a, maximum: $b")
# want one diagram to be alpha vs beta vs infected fraction 

# need to find transition or use a granular enough parameter regime