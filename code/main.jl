### GENERATE INPUT DATA
include("data-generation.jl")
#generate two sets of hypergraphs on 5k and 50k nodes
generate_data(5000,2,15)
generate_data(50000,2,15)

### GENERATE SIRS DATA 
#can add workers here or from the command line when initiating the environment.
# addprocs(25)
include("parallel-epidemics.jl")

### PROJECTED GRAPH DEGREES
include("degree-projection.jl")

#build figures
# include("figure1.jl")
# include("figure2.jl")
# include("figure3.jl")
# include("figure4_5.jl")
# include("figure6.jl")
#figure 7 
# include("spatial-graph-projected-lam1.jl")

