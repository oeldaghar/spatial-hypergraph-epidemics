#generate input data
include("data-generation.jl")
#starter data
generate_data(5000,2,15)
generate_data(50000,2,15)

### GENERATE DATA 
include("parallel-epidemics.jl")

#build figures.. manual process so far 
