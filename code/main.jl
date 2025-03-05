### GENERATE INPUT DATA
include("data-generation.jl")
#generate two sets of hypergraphs on 5k and 50k nodes
generate_data(5000,2,15)
generate_data(5000,5,15)
generate_data(50000,2,15)
generate_data(50000,5,15)

### GENERATE SIRS DATA 
#can add workers here or from the command line when initiating the environment.
# addprocs(25)
# include("epidemics/epidemic-experiments.jl")

# Generate additional data required for plotting specific figures.
# Data files are named as {file_path}/figure{num}-{figure_name}-data.jl
# Corresponding figure files are named {file_path}/figure{num}-{figure_name}.jl
include("code/figure5-alpha-figure-data.jl")
include("code/figure6-hypergraph_stats-data.jl")
include("code/clustering/figure7-hgcc-data.jl")
include("code/ppr/figure9-higher_order_ppr-data.jl")
include("code/figure15-lambda1-data.jl")
include("code/figure16-degree-data.jl")

# create figure path 
FIGURE_PATH = "data/output/figures/final/"
if !ispath(FIGURE_PATH)
    mkpath(FIGURE_PATH)
end

### generate figures 
include("code/figure1-hyperedge-formation.jl")
include("code/figure2-alpha-interpolation.jl")
include("code/figure3-distance-increasing.jl")
include("code/figure4-example-hypergraphs.jl")
# figures that require additional data. file names 'figure*-data.jl'
include("code/figure5-alpha-figure.jl")
include("code/figure6-hypergraph_stats.jl")
include("code/clustering/figure7-hgcc.jl")
include("code/ppr/figure9-higher_order_ppr.jl")
include("code/figure15-lambda1.jl")
include("code/figure16-degree.jl")
# back to figures where data already cached or not needed
include("code/clustering/figure8-cc-example.jl")
# epidemic figures. required aggregated data to be loaded first 
# load aggregated data 
include("code/epidemics/epidemic-figures-utils.jl")
# back to epidemic figures 
include("code/epidemics/figure10-total_infections_demo.jl")
include("code/epidemics/figure11-no_bistability.jl")
include("code/epidemics/figure12-total-infections.jl")
include("code/epidemics/figure13-hyperedge-effects-p1.jl")
include("code/epidemics/figure14-hyperedge-effects-p2.jl")
include("code/epidemics/figure17-hyperedge-transmissions.jl")