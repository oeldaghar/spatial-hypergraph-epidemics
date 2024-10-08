#generate input data
include("data-generation.jl")
#generate two sets of hypergraphs 
generate_data(5000,2,15)
generate_data(50000,2,15)

### GENERATE SIRS DATA 
# addprocs(25)
include("parallel-epidemics.jl")

#build figures
# include("figure1.jl")
# include("figure2.jl")
# include("figure3.jl")
#figures 4 + 5
# include("analysis.jl")
#figure 6 
# include("degree-projection.jl")
#figure 7 
# include("spatial-graph-projected-lam1.jl")

