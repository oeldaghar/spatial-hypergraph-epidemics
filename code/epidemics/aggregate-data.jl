using StatsBase
using JSON 
using JSON3
using ProgressMeter 

println("LOADING IN TOOLS...")
# screen 
include("code/data-io.jl")
include("code/hypergraph-tools.jl")
# REPL 
# include("../data-io.jl")
# include("../hypergraph-tools.jl")

data_dir = "data/epidemics/sirs/"
save_path = "aggregated_data/"
if !ispath(save_path)
    mkpath(save_path)
end

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

function load_epidemic_data(data_dir,fname;excluded=[],truncate=false,ntruncate=1000)
    excluded = Set(excluded)
    #initialize result 
    data = Dict{String,Array}()
    data["beta"] = Vector{Float64}()
    data["gamma"] = Vector{Float64}()
    data["delta"] = Vector{Float64}()
    data["exo"] = Vector{Float64}()
    data["tmax"] = Vector{Int}()
    data["hyper_beta_func"] = Vector{String}()
    data["initial_num_infected_nodes"] = Vector{Int}()
    data["infections"] = Vector{Vector{Float64}}()
    data["ninfected_hyperedge"] = Vector{Vector{Dict}}()
    data["ntransmissions_hyperedge"] = Vector{Vector{Dict}}()
    
    fpath = joinpath(data_dir,fname)
    println("LOADING $fname FROM $data_dir")
    open(fpath, "r") do file
        for (linecount,line) in enumerate(eachline(file))
            if linecount%100==0
                println("LOADED INFO FOR $linecount EPIDEMICS")
            end
            try
                obj = JSON3.read(line)
                for k in keys(obj)
                    key = string(k)
                    if !(key in excluded)
                        if truncate && typeof(obj[key])<:JSON3.Array
                            push!(data[key],obj[key][end-ntruncate:end])
                        else
                            push!(data[key],obj[key])
                        end
                    end
                end
            catch e
                # skip that row due to an error in parsing JSON
                println("Skipping row $linecount due to error: $e")
            end
        end
    end
    # post processing
    data["infections"] = reduce(hcat,data["infections"])'
    return data 
end

# similar to above but rather than load the entire data in, we reduce on the memory intensive keys 
# "ntransmissions_hyperedge" and "ninfected_hyperedge".
function load_epidemic_data_reduce(data_dir,fname,reduce_func;excluded=[],truncate=false,ntruncate=1000)
    excluded = Set(excluded)
    #initialize result 
    data = Dict{String,Array}()
    data["beta"] = Vector{Float64}()
    data["gamma"] = Vector{Float64}()
    data["delta"] = Vector{Float64}()
    data["exo"] = Vector{Float64}()
    data["tmax"] = Vector{Int}()
    data["hyper_beta_func"] = Vector{String}()
    data["initial_num_infected_nodes"] = Vector{Int}()
    data["infections"] = Vector{Vector{Float64}}()
    # keys on which to reduce during streaming. these are the bulk of the memory
    # data is originally Vector{Vector{Dict}} but reduce along one dim to make Vector{Dict}
    data["ninfected_hyperedge"] = Vector{Dict}()
    data["ntransmissions_hyperedge"] = Vector{Dict}()
    reduce_keys = Set(["ninfected_hyperedge","ntransmissions_hyperedge"])
    
    fpath = joinpath(data_dir,fname)
    println("LOADING $fname FROM $data_dir")
    open(fpath, "r") do file
        for (linecount,line) in enumerate(eachline(file))
            if linecount%100==0
                println("LOADED INFO FOR $linecount EPIDEMICS")
            end
            try 
                obj = JSON3.read(line)
                for k in keys(obj)
                    key = string(k)
                    if !(key in excluded || key in reduce_keys)
                        if truncate && typeof(obj[key])<:JSON3.Array 
                            push!(data[key],obj[key][end-ntruncate:end])
                        else
                            push!(data[key],obj[key])
                        end
                    elseif key in reduce_keys
                        # load it in, reduce and push reduced output to 
                        push!(data[key],reduce_func(obj[key]))
                    end
                end
            catch e
                println("Skipping row $linecount due to error: $e")
            end
        end
    end
    # post processing
    data["infections"] = reduce(hcat,data["infections"])'
    return data 
end

function average_over_dicts(vec_of_dicts::Vector{Dict})
    # take sum over all dicts 
    s_dict = Dict{Any,Float64}()
    for dict in vec_of_dicts
        for key in keys(dict)
            s_dict[key] = get(s_dict,key,0)+dict[key]
        end
    end
    nitems = lastindex(vec_of_dicts)
    for key in keys(s_dict)
        s_dict[key] = s_dict[key]/nitems
    end
    return s_dict
end

# not really a traililng average... technically this truncates and takes an average.
# should just add a truncation parameter to 'average_over_dicts' function above
function dictionary_trailing_average(vec_of_dicts::S, window::Int)::Dict where S
    # build a dict to track sum prior to average 
    trailing_avgs = Dict{Any, Float64}()
    for dict in vec_of_dicts[end-window:end]
        for (key,val) in pairs(dict)
            trailing_avgs[key] = get(trailing_avgs,key,0)+val
        end
    end
    
    for key in keys(trailing_avgs)
        trailing_avgs[key] = trailing_avgs[key]/window
    end
    return trailing_avgs
end

function rolling_hyperedge_stats(data,hyperedge_key="ntransmissions_hyperedge",stream=false)
    # group data by starting condition, beta, and hyper_beta_func
    filtered_data = hcat(data["beta"], data["hyper_beta_func"])
    grouped_inds = group_by_function(x -> (x[1],x[2]), eachrow(filtered_data))
    # compute rolling averages for each trajectory/simulation and average over group size 
    println("Computing averages...")
    newdata = Dict()
    for (key,inds) in pairs(grouped_inds)
        # for original data 
        if !stream
            rolling_avgs = map(ind->dictionary_trailing_average(data[hyperedge_key][ind],1000),inds)
            group_rolling_avg = average_over_dicts(rolling_avgs)
        else
            group_rolling_avg = average_over_dicts(data[hyperedge_key][inds])
        end
        newdata[key] = group_rolling_avg
    end

    # get max hyperedge size 
    hyperedge_keys = string.(keys(newdata[first(keys(newdata))]))
    hmin,hmax = extrema(parse.(Int,hyperedge_keys))
    
    beta_vals = sort(unique(data["beta"]))
    pdata = Dict()
    for hyper_beta_func in unique(data["hyper_beta_func"])
        # main plotting setup 
        plotting_data = zeros(Float64,hmax,lastindex(beta_vals))
        for (i,beta) in enumerate(beta_vals)
            dict = newdata[(beta,hyper_beta_func)] # filter to hyper_beta_func 
            for key in keys(dict)
                plotting_data[parse(Int,string(key)),i] = dict[key]
            end
        end

        pdata[hyper_beta_func] = plotting_data
    end
    return pdata 
end

function hypergraph_stats_data(data_dir,fname,stream=false)
    # load data 
    if stream 
        data = load_epidemic_data_reduce(data_dir,fname,x->dictionary_trailing_average(x,1000),truncate=true,ntruncate=1000)
    else
        data = load_epidemic_data(data_dir,fname)    
    end
    # generate individual figures 
    pdata = Dict()
    for hyperkey in ["ntransmissions_hyperedge","ninfected_hyperedge"]
        output = rolling_hyperedge_stats(data,hyperkey,stream)
        pdata[hyperkey] = output 
    end
    # handle infections data 
    filtered_data = hcat(data["beta"], data["hyper_beta_func"])
    grouped_inds = group_by_function(x -> (x[1], x[2]), eachrow(filtered_data))
    rolling_avg_infections = Dict(k => rolling_average(data["infections"][grouped_inds[k], :], 1000) for k in keys(grouped_inds))

    static_keys = collect(keys(rolling_avg_infections))
    regrouped_keys = group_by_function(x->x[2], static_keys)
    infs_dict = Dict()
    for key in keys(regrouped_keys)
        sorted_betas = sort(static_keys[regrouped_keys[key]])
        infs_dict[key] = [rolling_avg_infections[x] for x in sorted_betas]
    end
     
    pdata["infections"] = infs_dict
    return pdata
end

# for padding matricies (for use with heapmaps of different sizes)
function _pad_matrix(matrix, nrows=100)
    nnewrows = nrows - size(matrix, 1)
    ncols = size(matrix, 2)
    if nnewrows > 0
        matrix = vcat(matrix, zeros(nnewrows, ncols))
    end
    return matrix
end

# fraction of transmissions vs hyperedge size (column normalization)
# infection efficiency: fraction of transmissions vs fraction hyperedges  (column normalization / hyperedge fraction)
# fraction hyperedges vs fraction coverage (hyperedge fraction * fraction of nodes covered)
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

beta_vals = vcat(1e-3:1e-3:1e-2, 2e-2:1e-2:1e-1, 2e-1:1e-1:9e-1)
beta_vals = sort(unique(beta_vals))

base_graph_names = [
    # "spatial-hypergraph-5000-2",
    # "spatial-hypergraph-5000-5",
    # larger 
    "spatial-hypergraph-50000-2", # graph we have all epidemics for 
    # "spatial-hypergraph-50000-5"
]

data_paths = [
    "data/epidemics/sirs/scratch-v1/"
    "data/epidemics/sirs/"
    "data/epidemics/sirs/scratch-v2/"
    "data/epidemics/sirs/scratch-v3/"
]
for data_path in data_paths
    # LOAD + AGGREGATE DATA
    RAW_DATA = Dict()
    for base_graph in base_graph_names 
        graph_names = map(x->"$(splitext(basename(x))[1]).jsonl",get_fnames(base_graph))
        filter!(x->occursin("newalpha.jsonl",x),graph_names)

        println("ALL GRAPHS..")
        for x in graph_names
            println(x)
        end
        for fname in graph_names 
            RAW_DATA[fname] = load_epidemic_data_reduce(data_path,fname,x->dictionary_trailing_average(x,1000),truncate=false)
        end
    end 
    # SAVE DATA 
    last_dir = split(data_path,"/")[end-1]
    fpath = ""
    if last_dir=="sirs"
        fpath = joinpath(save_path,"aggregated-sirs-output.json")
    else
        fpath = joinpath(save_path,"aggregated-sirs-output-$last_dir.json")
    end
    open(fpath,"w") do io
        JSON3.write(io, RAW_DATA)
    end     
end
