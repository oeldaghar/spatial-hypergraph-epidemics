using StatsPlots
using Plots
using CairoMakie
using Measures
using StatsBase
using JSON 
using JSON3
using ProgressMeter
import LinearAlgebra: norm
using LaTeXStrings
using DelimitedFiles

println("LOADING IN TOOLS...")
include("../data-io.jl")
include("../hypergraph-tools.jl")

data_dir = "data/epidemics/sirs/"

## GLOBALS 
ALPHA_VALS = collect(range(0,2,15))
BETA_VALS = vcat(1e-3:1e-3:1e-2, 2e-2:1e-2:1e-1, 2e-1:1e-1:9e-1)
BETA_VALS = sort(unique(BETA_VALS))

IND_TO_BETA = Dict(enumerate(BETA_VALS))
BETA_TO_IND = Dict(v=>k for (k,v) in pairs(IND_TO_BETA))

# helper functions 
function rolling_average(input,window)
    return mean(input[end-window:end])
end

# values of keys are inds with same function evaluation
function group_by_function(f, itr)
    result = Dict()
    for (i, v) in enumerate(itr)
        key = f(v)
        push!(get!(result, key, Int[]), i)
    end
    return result
end

# converts a dict into a vector. useful for parsing the ntransmissions_hyperedge and ninfected_hyperedge vector of dicts
function make_matrix(mydict::Dict, max_key = -1)
    if max_key<0
        max_key = maximum(parse.(Int, keys(mydict)))
    end
    result = zeros(max_key)
    for (key, value) in pairs(mydict)
        result[parse(Int, key)] = value
    end
    return result
end

function hyperedge_information(fname::String)
    nnodes = parse(Float64,split(fname,"-")[3])
    hedges = read_hedges(fname)
    # edge coverage by size 
    hyper_sizes = map(x->length(x),hedges)
    max_hyper_size = maximum(hyper_sizes)
    edge_coverage = zeros(max_hyper_size)
    for val in hyper_sizes
        edge_coverage[val] +=1
    end
    # node coverage by size # node_coverage[k] = union(hyperedges of size k) 
    binned_hyper_edges = Dict{Int,Vector{Int}}()
    for (ind,val) in enumerate(hyper_sizes)
        if haskey(binned_hyper_edges,val)
            push!(binned_hyper_edges[val],ind)
        else
            binned_hyper_edges[val] = Vector{Int}([ind])
        end
    end
    node_coverage = zeros(max_hyper_size)
    for (hyper_size,inds) in pairs(binned_hyper_edges)
        node_coverage[hyper_size] = length(reduce(union,hedges[inds]))
    end
    return node_coverage, edge_coverage
end

function load_aggregated_data(data_dir,fname)
    start_time = time()
    raw_data = open(joinpath(data_dir,fname), "r") do io
        JSON3.read(io, Dict)
    end
    end_time = time()
    println("loaded data in $(end_time-start_time) seconds")
    println("reshaping data...")
    start_time = time()
    for key in keys(raw_data)
        nepidemics = lastindex(raw_data[key]["beta"])

        num_trailing_steps = lastindex(raw_data[key]["infections"])/nepidemics
        num_trailing_steps = Int(num_trailing_steps)

        raw_data[key]["infections"] = reshape(raw_data[key]["infections"],(nepidemics,num_trailing_steps))
    end
    for key in keys(raw_data)
        max_hsize = maximum(map(x->maximum(parse.(Int,keys(x))),raw_data[key]["ninfected_hyperedge"]))
        raw_data[key]["ninfected_hyperedge"] = Array(reduce(hcat, make_matrix.(raw_data[key]["ninfected_hyperedge"],max_hsize))')
        raw_data[key]["ntransmissions_hyperedge"] = Array(reduce(hcat, make_matrix.(raw_data[key]["ntransmissions_hyperedge"],max_hsize))')
        # correct over counts in ninfected_hyperedge 
        h_correction = collect(1:max_hsize).-1
        h_correction[1] = 1 # should be 0 but numerator for this is 0 so result is okay
        
        raw_data[key]["ninfected_hyperedge"] = raw_data[key]["ninfected_hyperedge"] ./ h_correction'
    end
    end_time = time()
    println("reshaped data in $(end_time-start_time) seconds")
    return raw_data
end
load_aggregated_data(fname) = load_aggregated_data("",fname)

# load in aggregated data.. takes about 5mins
aggregated_data = load_aggregated_data("aggregated_data","aggregated-sirs-output.json")
aggregated_data1 = load_aggregated_data("aggregated_data","aggregated-sirs-output-scratch-v1.json")
aggregated_data2 = load_aggregated_data("aggregated_data","aggregated-sirs-output-scratch-v2.json")
aggregated_data3 = load_aggregated_data("aggregated_data","aggregated-sirs-output-scratch-v3.json")
