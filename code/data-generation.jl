#code for generating data. 
#will generate a spatial hypergraph and then it's pairwise projection
using Random, MatrixNetworks
include("hypergraph-tools.jl")
include("spatial-hypergraph.jl")
using DelimitedFiles
using ProgressMeter

"""
    generate_data(n,k;rseed=1,ngraphs=5)

generates k hypergraphs on n nodes as defined by spatial_hypergraph_edges(n,k,degreedist)
"""
function generate_data(n=5000,d=2,ngraphs=5;rseed=1)
    #interpolates between the alpha=0 (pairwise case) and alpha=2 (entire neighborhood as hyperedge)
    alphas = range(0,2,length=ngraphs)

    println("GENERATING SYNTHETIC GRAPHS..")
    @showprogress for i=1:ngraphs
        alpha = alphas[i]
        #main generation call 
        Random.seed!(rseed)
        X = rand(d,n) #d = dimension, n = num nodes  
        degreedist=LogNormal(log(3),1)
        degs = rand(degreedist,n)
        degs = min.(ceil.(Int,degs),n-1)
        hedges,xy = hypergraph_edges(X;degs=degs,radfunc=get_func(alpha))
        #cleaning up 
        sort!.(hedges)
        sort!(hedges)
        unique!(hedges)
        #save info 
        gname = "spatial-hypergraph-$n-$d-$alpha-$rseed-newalpha"
        file_name = "data/hypergraphs/$gname.txt"
        save_hedges(hedges,file_name)
        file_name = "data/hypergraphs/$gname.xy"
        open(file_name, "w") do io
            writedlm(io, xy')
        end
    end
end