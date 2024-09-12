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
function generate_data(n=5000,d=2,ngraphs=5;rseed=1)
    Random.seed!(rseed)
    random_seeds = abs.(rand(Int,ngraphs))

    #interpolates between the alpha=0 (pairwise case) and alpha=2 (entire neighborhood as hyperedge)
    alphas = range(0,2,length=ngraphs)

    for i=1:ngraphs
        random_seed = random_seeds[i]
        alpha = alphas[i]
        #set random seed 
        Random.seed!(random_seed)
        #main generation call 
        X = rand(d,n) #d = dimension, n = num nodes  
        hedges,xy = hypergraph_edges(X;degreedist=LogNormal(log(3),1),radfunc=get_func(alpha))
        #cleaning up 
        sort!.(hedges)
        sort!(hedges)
        unique!(hedges)
        #save info 
        file_name = "data/hypergraphs/spatial-hypergraph-$n-$d-$alpha-$random_seed.txt"
        save_hedges(hedges,file_name)
        file_name = "data/hypergraphs/spatial-hypergraph-$n-$d-$alpha-$random_seed.xy"
        open(file_name, "w") do io
            writedlm(io, xy')
        end
    end
end