using Distributed
addprocs(25)

@everywhere include("hypergraph-tools.jl")
@everywhere include("spatial-hypergraph.jl")
