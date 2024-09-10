#code for generating data. 
#will generate a spatial hypergraph and then it's pairwise projection

using Random, MatrixNetworks
include("hypergraph-tools.jl")
include("spatial-hypergraph.jl")
using DelimitedFiles

"""
    generate_data(n,k;rseed=1,ngraphs=5)

generates k hypergraphs on n nodes as defined by spatial_hypergraph_edges(n,k,degreedist)
"""
function generate_data(n,k,ngraphs=5;rseed=1)
    Random.seed!(rseed)
    random_seeds = abs.(rand(Int,ngraphs))

    for random_seed in random_seeds
        Random.seed!(random_seed)
        hedges,xy = spatial_hypergraph_edges(n,k,degreedist=LogNormal(log(3),1))
        #preprocess
        sort!.(hedges)
        sort!(hedges)
        unique!(hedges)

        #save info 
        file_name = "data/hypergraphs/spatial-hypergraph-$n-$random_seed.txt"
        save_hedges(hedges,file_name)

        file_name = "data/hypergraphs/spatial-hypergraph-$n-$random_seed.xy"
        open(file_name, "w") do io
            writedlm(io, xy')
        end
    end
end