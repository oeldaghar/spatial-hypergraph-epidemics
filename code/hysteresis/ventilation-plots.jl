# initial figures for hysteresis effects 
using StatsPlots
using Plots
using Measures
using StatsBase
using JSON 
using JSON3
using ProgressMeter

println("LOADING IN TOOLS...")
include("../data-io.jl")
include("../hypergraph-tools.jl")

data_dir = "data/hysteresis/sirs/"
figure_path = "data/hysteresis/figures"
if !ispath(figure_path)
    mkpath(figure_path)
end

original_graph_names = [
            # 5k nodes in 2d
            ["spatial-hypergraph-5000-2-0.0-1-newalpha.jsonl",
                "spatial-hypergraph-5000-2-1.0-1-newalpha.jsonl",
                "spatial-hypergraph-5000-2-2.0-1-newalpha.jsonl"],
            # 5k nodes in 5d
            ["spatial-hypergraph-5000-5-0.0-1-newalpha.jsonl",
                "spatial-hypergraph-5000-5-1.0-1-newalpha.jsonl",
                "spatial-hypergraph-5000-5-2.0-1-newalpha.jsonl"],
            # 50k nodes in 2d
            ["spatial-hypergraph-50000-2-0.0-1-newalpha.jsonl",
                "spatial-hypergraph-50000-2-1.0-1-newalpha.jsonl",
                "spatial-hypergraph-50000-2-2.0-1-newalpha.jsonl"],
            # 50k nodes in 5d
            ["spatial-hypergraph-50000-5-0.0-1-newalpha.jsonl",
                "spatial-hypergraph-50000-5-1.0-1-newalpha.jsonl",
                "spatial-hypergraph-50000-5-2.0-1-newalpha.jsonl"],
]

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

# data_dir = "data/hysteresis/sirs/"
# fname = "spatial-hypergraph-5000-5-2.0-1-newalpha.jsonl"
# fpath = joinpath(data_dir,fname)
# data = load_epidemic_data(data_dir,fname,excluded=["ntransmissions_hyperedge","ninfected_hyperedge"],truncate=true)
# data["ntransmissions_hyperedge"] #should return empty 

# plot 1: can we see hysteresis as we vary beta for a single graph 
# plot binned rolling average by starting condition vs beta as scatter plot
# if there is no bistability (not in that parameter regime), then
# we'd expect to see that points corresponding to the same beta should be very close to one another
# (the initial condition doesn't matter)
# otherwise, we'd expect to see meaningfully different values for the same beta 

# Example usage
function plot_hysteresis(data_dir, fname)
    data = load_epidemic_data(data_dir,fname,
                                excluded = ["ntransmissions_hyperedge","ninfected_hyperedge"],
                                truncate=true,
                                ntruncate=1000,
                            )

    # group data by starting condition, beta, and hyper_beta_func
    filtered_data = hcat(data["beta"], data["hyper_beta_func"], data["initial_num_infected_nodes"])
    grouped_inds = group_by_function(x -> (x[1], x[2], x[3]), eachrow(filtered_data))
    grouped_rolling_averages = Dict(k => rolling_average(data["infections"][grouped_inds[k], :], 1000) for k in keys(grouped_inds))

    # put the plot together 
    static_keys = sort(collect(keys(grouped_rolling_averages)))
    linear_keys = [k for k in static_keys if k[2] == "linear"]
    sqrt_keys = [k for k in static_keys if k[2] == "sqrt"]
    pairwise_keys = [k for k in static_keys if k[2] == "pairwise"]

    beta_vals = unique(first.(linear_keys))
    beta_map = Dict(map(x -> (x[2], x[1]), enumerate(sort(beta_vals))))

    # beta_map[beta_val] = x_ind
    function _plot_keys(keys, title_suffix)
        # xs = [beta_map[x[1]] for x in keys]
        xs = [x[1] for x in keys]
        ys = [grouped_rolling_averages[x] for x in keys]

        f = Plots.scatter(xs, ys, leg = false,
                xticks = ([1e-3,1e-2,1e-1,1e0], ["1e-3","1e-2","1e-1","1e0"]),
                xlabel = "Pairwise Infection Probability", 
                ylabel = "Binned Average\nRolling Infections", 
                title = "Beta vs Rolling Infections - $title_suffix\n$fname",
                markerstrokewidth = 0,
                xlims=(8e-4,1.05e0))
        Plots.plot!(f,xscale=:log10)
    end

    return _plot_keys(linear_keys, "Linear"),_plot_keys(sqrt_keys, "Sqrt"), _plot_keys(pairwise_keys, "Pairwise")
end

################################
## Figure 1 - hysteresis testing
################################
# Example usage
# data_dir = "data/hysteresis/sirs/"
# fname = "spatial-hypergraph-5000-2-2.0-1-newalpha.jsonl"
# f1,f2,f3 = plot_hysteresis(data_dir, fname)
# display(f1)
# display(f2)
# display(f3)

### 
# fs = []
# for fname in [
#             "spatial-hypergraph-5000-2-0.0-1-counter_1.json",
#             "spatial-hypergraph-5000-2-1.0-1-counter_1.json",
#             "spatial-hypergraph-5000-2-2.0-1-counter_1.json"
#         ]
#     println("WORKING ON $fname")
#     f1,f2 = plot_hysteresis(data_dir, fname)
#     ylims!(f1,(-50,3000))
#     ylims!(f2,(-50,3000))
#     push!(fs,f1)
#     push!(fs,f2)
# end

# # side by side plot as we vary alpha across columns
# f_main = plot(fs[1], fs[3], fs[5], fs[2], fs[4], fs[6],
#               layout=(2, 3),
#               size=(1000, 600),
#               left_margin=5mm, right_margin=2mm, top_margin=2mm, bottom_margin=5mm,
#               title=["alpha=0" "alpha=1" "alpha=2" "alpha=0" "alpha=1" "alpha=2"],
#               titlefontsize=12,
#               plot_title="Linear Normalization",
#               plot_titlefontsize=18)
# # remove xlabel, xtick labels, and adjust spacing for first row 
# for i=1:3
#     plot!(f_main[i], xticks=(xticks(f_main[i])[1], ""),
#         xlabel="",
#         bottom_margin=8mm,
#         top_margin=8mm)
# end
# # remove ylabel for 2nd column and beyond
# for i=1:6
#     if i%3!=1
#         plot!(f_main[i], yticks=(yticks(f_main[i])[1], ""),
#             ylabel="",
#             left_margin=3mm)
#     end
# end
# # Add row titles using annotation
# plot!(f_main, annotation=(5e-2, 3800, text("Square Root Normalization", 18, :center)),
#                 subplot=5)
# # save figure 
# savefig(f_main, joinpath(figure_path,"hysteresis_plot_v1.pdf"))
# plot!(f_main,dpi=1000)
# savefig(f_main, joinpath(figure_path,"hysteresis_plot_v1.png"))

#zooming in on one piece 
function stochastic_evolution(data_dir,fname)
    data = load_epidemic_data(data_dir,fname,excluded=["ntransmissions_hyperedge","ninfected_hyperedge"])

    # group data by starting condition, beta, and hyper_beta_func
    filtered_data = hcat(data["beta"], data["hyper_beta_func"], data["infections"])
    grouped_inds = group_by_function(x -> (x[1], x[2], x[3]), eachrow(filtered_data))

    # put the plot together 
    static_keys = sort(collect(keys(grouped_inds)))
    linear_keys = [k for k in static_keys if k[2] == "linear"]
    sqrt_keys = [k for k in static_keys if k[2] == "sqrt"]

    beta_vals = unique(first.(linear_keys))
    beta_map = Dict(map(x -> (x[2], x[1]), enumerate(sort(beta_vals))))

    beta = 9e-3
    # pick a value of beta and do a plot for each normalization
    inds = (first.(static_keys).==beta) .& ([x[2]=="linear" for x in static_keys])

    f = Plots.plot()
    for (i,key) in enumerate(static_keys[inds])
        group = grouped_inds[key]
        Plots.plot!(f,data["infections"][group,1:200]',c=i,leg=false)
    end

    alpha_val = split(fname,"-")[5]
    Plots.plot!(f,
        xlabel="Time Step",
        ylabel="Total Infections",
        title = "alpha=$alpha_val, beta = $beta\n$fname")
    return f 
end

# f1 = stochastic_evolution(data_dir,fname)

# filtered_fnames = [
#             "spatial-hypergraph-5000-2-0.0-1-counter_1.json",
#             "spatial-hypergraph-5000-2-1.0-1-counter_1.json",
#             "spatial-hypergraph-5000-2-2.0-1-counter_1.json"
#         ]
# fs = []
# for fname in filtered_fnames
#     f = stochastic_evolution(data_dir,fname)
#     ylims!(f,(0,5000))
#     push!(fs,f)
# end

# #put them together 
# main_f = plot(fs...,layout=(3,1),
#         size = (250,650),
#         # margins=8mm,
#         left_margin = 7mm,
#         titlefontsize=11)
# savefig(main_f, joinpath(figure_path,"hysteresis_stochastic_plot_v1.pdf"))
# plot!(main_f,dpi=1000)
# savefig(main_f, joinpath(figure_path,"hysteresis_stochastic_plot_v1.png"))
        
#############################
## PLOTS for hypergraph stats 
#############################
# average over dicts 
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
    # handle infections 
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

# plot individual figures 
function make_transmission_heatmaps(data_dict,hyperedge_dict,column_normalization=false)
    fig_dict = Dict()
    for (key,plotting_data) in pairs(data_dict)
        alpha,hyperkey,hyper_beta_func = key
        curr_plotting_data = deepcopy(plotting_data)

        beta_tick_map = Dict(reverse.(enumerate(beta_vals)))
        ticks = (vcat([beta_tick_map[x] for x in [1e-3,1e-2,1e-1]],lastindex(beta_vals)) ,["1e-3","1e-2","1e-1","1e0"])
        if hyperkey == "ntransmissions_hyperedge"
            mylabel = "Rolling Average Transmission\n by Hyperedge Size"
        elseif hyperkey == "ninfected_hyperedge"
            mylabel = "Rolling Average Infection\n by Hyperedge Size"
            # recorded sum_v sum_{h: v\in h} (number of infected nodes adjacent to v) 
            # this overcounts by |h|-1 for each hyperedge (split into inf and not inf and count)
            curr_plotting_data = mapslices(x->x./(0:lastindex(x)-1),curr_plotting_data,dims=1)
        end
        
        hinfo = deepcopy(hyperedge_dict[alpha]["node_coverage"])
        nnodes = hyperedge_dict[0.0]["node_coverage"][2,2] # node coverage by pairwise in pairwise 
        hinfo = mapslices(x->x./nnodes,hinfo,dims=1)
        curr_plotting_data ./= hinfo

        if column_normalization
            curr_plotting_data = mapslices(x -> x ./ maximum(x[(!).(isnan.(x))]), curr_plotting_data, dims=1)
            # removes portions of the heatmap with no data
            curr_plotting_data = (10.0).^log10.(curr_plotting_data)
        else
            curr_plotting_data = 1 .+log10.(curr_plotting_data)
        end
        
        f = Plots.heatmap(curr_plotting_data,
                yscale=:log10,            
                cmap=:viridis,
                ylims=(0.9,size(curr_plotting_data,1)+10),
                xticks = ticks,
                xlabel = "Pairwise Infection Probability",
                ylabel = mylabel,
                title="Alpha:$alpha Norm: $hyper_beta_func")
        fig_dict[key] = f
    end
    return fig_dict
end

beta_vals = vcat(1e-3:1e-3:1e-2, 2e-2:1e-2:1e-1, 2e-1:1e-1:9e-1)
beta_vals = sort(unique(beta_vals))

base_graph_names = [
    "spatial-hypergraph-5000-2",
    "spatial-hypergraph-5000-5",
    # larger 
    "spatial-hypergraph-50000-2",
    "spatial-hypergraph-50000-5"
]

function get_hypergraph_info(fnames,max_size)
    hypergraph_info = Dict()
    for fname in fnames     
        node_coverage,edge_coverage = hyperedge_information(fname)
        node_coverage = repeat(node_coverage,outer = (1,lastindex(beta_vals)))
        node_coverage = _pad_matrix(node_coverage,max_size)

        edge_coverage = repeat(edge_coverage,outer = (1,lastindex(beta_vals)))
        edge_coverage = _pad_matrix(edge_coverage,max_size)
        alph = _get_alpha(fname)
        gdict = Dict()
        gdict["node_coverage"] = node_coverage
        gdict["edge_coverage"] = edge_coverage
        hypergraph_info[alph] = gdict
    end
    return hypergraph_info 
end

base_graph_names
gnames = get_fnames(base_graph_names[1],"newalpha")

# load cached data dict 
data_dir = "data/hysteresis/figures/"
fnames = [
    "spatial-hypergraph-50000-2-aggregated_data.json",
    "spatial-hypergraph-50000-5-aggregated_data.json",
    "spatial-hypergraph-5000-2-aggregated_data.json",
    "spatial-hypergraph-5000-5-aggregated_data.json",
]

fname = fnames[3]
base_graph = String(first(split(fname,"-aggregated_data")))
graphs = get_fnames(base_graph,"newalpha")

# Read the JSON file as a string
json_string = read(joinpath(data_dir,fname), String)
# Parse the JSON string
agg_inf_data = JSON3.read(json_string)

function parse_key(key::Symbol)
    # Convert the symbol to a string and remove the outer parentheses
    str = String(key)[2:end-1]
    
    # Split the string by comma and space
    parts = split(str, ", ")
    
    # Parse the float
    float_val = parse(Float64, parts[1])
    
    # Remove the quotes from the strings
    string1 = parts[2][2:end-1]
    string2 = parts[3][2:end-1]
    
    # Return the tuple
    return (float_val, String(string1), String(string2))
end

# Example usage:
# key = Symbol("(0.42857142857142855, \"ninfected_hyperedge\", \"linear\")")
# result = parse_key(key)
# println(result)  # Output: (0.42857142857142855, "ninfected_hyperedge", "linear")
# println(typeof(result))  # Output: Tuple{Float64, String, String}
agg_inf_data = Dict(parse_key(k)=>Array(agg_inf_data[k]) for k in keys(agg_inf_data))
hyperedge_keys = [x for x in keys(agg_inf_data) if occursin("hyperedge",x[2])]
# reshape infection data 
max_hyperedge_size = 0
for key in hyperedge_keys    
    max_hyperedge_size = max(max_hyperedge_size,Int(size(agg_inf_data[key],1)/length(beta_vals)))
end
mat_shape = (max_hyperedge_size,length(beta_vals))

h_inf_data = Dict{Tuple,Array}()
for key in keys(agg_inf_data)
    if key in hyperedge_keys
        h_inf_data[key] = reshape(agg_inf_data[key],mat_shape)
        # remove key for inf data 
        delete!(agg_inf_data,key)
    end
end

# compute hyperedge information 
hypergraph_info = get_hypergraph_info(graphs,max_hyperedge_size)
# rebuild heatmaps (as test)
# figs = make_transmission_heatmaps(agg_inf_data,hypergraph_info,true)
# figs[2.0,"ntransmissions_hyperedge","pairwise"]
# plot the relationship we care about 

# bin the data into coarser bins 
function _custom_heatmap(pdata)
    # function to bin rows of the heatmap weight matrix
    function _bin_indices(bins::Vector, pdata::Vector)
        bin_dict = Dict{Float64, Vector{Int}}()
        for i in 1:length(bins)-1
            bin_dict[bins[i]] = findall(x -> bins[i] <= x < bins[i+1], pdata)
        end
        bin_dict[bins[end]] = findall(x -> bins[end] <= x, pdata)
        return bin_dict
    end
    
    function _log_bin_ydata(pdata,ybins)
        binned_inds = _bin_indices(ybins,collect(1:lastindex(pdata,1)))
        new_mat = zeros(Float64,lastindex(ybins)-1,lastindex(pdata,2))
        for (newrow,key) in enumerate(ybins[1:end-1])
            old_rows = binned_inds[key]
            # account for nan values occuring in different rows/cols in sum 
            for col=1:lastindex(pdata,2)
                tmp = @view pdata[old_rows,col]
                # get non-nan values and sum 
                inds = (!).(isnan.(tmp))
                if length(inds)>0
                    new_mat[newrow,col] = sum(tmp[inds])
                end
                # new_mat[newrow,:] = sum(x->!isnan(x),pdata[old_rows,:],dims=1)
            end
        end
        return new_mat
    end
    
    max_hsize = lastindex(pdata,1)
    ybins = (10.0).^(range(1,log10(max_hsize+10),21))
    ybins = vcat(1:9,ybins)
    ybins = sort(unique(ybins))

    new_mat = _log_bin_ydata(pdata,ybins)
    # renormalize columns
    new_mat = mapslices(x -> x ./ maximum(x[(!).(isnan.(x))]), new_mat, dims=1)
    # yvals = [mean(ybins[i:i+1]) for i=1:lastindex(ybins)-1]
    # xrange, yrange, data_matrix, yaxis scale
    f = Plots.heatmap(1:lastindex(new_mat,2), ybins[1:end-1], log10.(new_mat),
    # f = Plots.heatmap(1:lastindex(new_mat,2), yvals, log10.(new_mat),
                yscale=:log10,
                color=:viridis,
                clims=(-5,0),
                )
    return f#,new_mat,ybins 
end

function new_transmission_heatmaps(data_dict,hyperedge_dict,column_normalization=false)
    fig_dict = Dict()
    for (key,plotting_data) in pairs(data_dict)
        alpha,hyperkey,hyper_beta_func = key
        curr_plotting_data = deepcopy(plotting_data)

        beta_tick_map = Dict(reverse.(enumerate(beta_vals)))
        ticks = (vcat([beta_tick_map[x] for x in [1e-3,1e-2,1e-1]],lastindex(beta_vals)) ,["1e-3","1e-2","1e-1","1e0"])
        if hyperkey == "ntransmissions_hyperedge"
            mylabel = "Rolling Average Transmission\n by Hyperedge Size"
        elseif hyperkey == "ninfected_hyperedge"
            mylabel = "Rolling Average Infection\n by Hyperedge Size"
            # recorded sum_v sum_{h: v\in h} (number of infected nodes adjacent to v) 
            # this overcounts by |h|-1 for each hyperedge (split into inf and not inf and count)
            curr_plotting_data = mapslices(x->x./(0:lastindex(x)-1),curr_plotting_data,dims=1)
        end
        
        hinfo = deepcopy(hyperedge_dict[alpha]["node_coverage"])
        nnodes = hyperedge_dict[0.0]["node_coverage"][2,2] # node coverage by pairwise in pairwise 
        hinfo = mapslices(x->x./nnodes,hinfo,dims=1)
        curr_plotting_data ./= hinfo

        if column_normalization
            curr_plotting_data = mapslices(x -> x ./ maximum(x[(!).(isnan.(x))]), curr_plotting_data, dims=1)
            # removes portions of the heatmap with no data
            curr_plotting_data = (10.0).^log10.(curr_plotting_data)
        # else
        #     curr_plotting_data = 1 .+log10.(curr_plotting_data)
        end
        
        f = _custom_heatmap(curr_plotting_data)
        Plots.plot!(f,
                cmap=:viridis,
                ylims=(0.9,size(curr_plotting_data,1)+10),
                xticks = ticks,
                xlabel = "Pairwise Infection Probability",
                ylabel = mylabel,
                title="Alpha:$alpha Norm: $hyper_beta_func"
        )
        fig_dict[key] = f
    end
    return fig_dict
end

# new_figs = new_transmission_heatmaps(agg_inf_data, hypergraph_info)
# new_figs[(2.0,"ntransmissions_hyperedge","sqrt")]

function node_normalized_stats(inf_dict, hypergraph_dict)
    new_inf_dict = deepcopy(inf_dict)
    for key in keys(new_inf_dict)
        alpha,hyperkey,hyper_beta_func = key
        curr_data = new_inf_dict[key]
        if hyperkey == "ninfected_hyperedge"
            # recorded sum_v sum_{h: v\in h} (number of infected nodes adjacent to v) 
            # this overcounts by |h|-1 for each hyperedge (split into inf and not inf and count)
            curr_data = mapslices(x->x./(0:lastindex(x)-1),curr_data,dims=1)
        end

        hinfo = deepcopy(hypergraph_dict[alpha]["node_coverage"])
        nnodes = hypergraph_dict[0.0]["node_coverage"][2,2] # node coverage by pairwise in pairwise 
        hinfo = mapslices(x->x./nnodes,hinfo,dims=1)
        new_inf_dict[key] ./= hinfo
    end
    return new_inf_dict
end

function nanfindmax(x::Vector)
    inds = (!).(isnan.(x))
    y = x[inds]
    maxval,maxind = findmax(y)
    new_inds = cumsum(inds) # index in the filtered data
    return findfirst(new_inds.==maxind)
end
nn_h_inf_data = node_normalized_stats(h_inf_data,hypergraph_info)

# tmp = deepcopy(nn_h_inf_data[(2.0,"ntransmissions_hyperedge","sqrt")])
# tmp = mapslices(x -> x ./ maximum(x[(!).(isnan.(x))]), tmp, dims=1)
# # removes portions of the heatmap with no data
# tmp = (10.0).^log10.(tmp)
# f = _custom_heatmap(tmp)


####### PLOTS #############
function make_transmission_plots(nn_h_inf_data,agg_inf_data)
    for hyper_beta_func in ["pairwise","linear","sqrt"]
        for hyper_stat in ["ntransmissions_hyperedge","ninfected_hyperedge"]

        end
    end
end
# hyper_beta_func = "pairwise"
# hyper_stat = "ntransmissions_hyperedge"

function beta_hyperstat_plot(nn_h_inf_data,hyper_stat,hyper_beta_func)
    f = Plots.plot(xlabel=L"\beta", ylabel="Hyperedge Size",
        title = "$hyper_stat\n$hyper_beta_func normalization")
    colors = cgrad(:viridis, 15, categorical=true)
    for (ind,alph) in enumerate(range(0,2,15))    
        tmp = deepcopy(nn_h_inf_data[(alph,hyper_stat,hyper_beta_func)])
        Plots.plot!(f, beta_vals,vec(mapslices(x -> nanfindmax(x), tmp, dims=1)),
                yscale=:log10,
                xscale=:log10,
                c=colors[ind], linewidth=1.2, leg = false, 
                # alpha=exp(alph)/exp(2)
                )
    end
    return f
end

using CairoMakie
using LaTeXStrings


beta_hyperstat_plot(nn_h_inf_data,"ntransmissions_hyperedge","linear")
beta_hyperstat_plot(nn_h_inf_data,"ntransmissions_hyperedge","sqrt")
beta_hyperstat_plot(nn_h_inf_data,"ntransmissions_hyperedge","pairwise")


function alpha_hyperstat_plot(nn_h_inf_data,hyper_stat,hyper_beta_func)
    f = Plots.plot(xlabel=L"\alpha", ylabel="Hyperedge Size",
        title = "$hyper_stat\n$hyper_beta_func normalization")
    for (ind,beta) in enumerate(beta_vals)    
        a_data = reduce(hcat,[nn_h_inf_data[(a,hyper_stat,hyper_beta_func)][:,ind] for a in range(0,2,15)])
        transparency_range = 4
        Plots.plot!(f, range(0,2,15),vec(mapslices(x -> nanfindmax(x), a_data, dims=1)),
                yscale=:log10,    
                c=3, linewidth=1.8, leg = false,
                alpha = exp(ind/transparency_range)/exp(length(beta_vals)/transparency_range)
            )         
    end
    f

end

alpha_hyperstat_plot(nn_h_inf_data,"ntransmissions_hyperedge","linear")
alpha_hyperstat_plot(nn_h_inf_data,"ntransmissions_hyperedge","sqrt")
alpha_hyperstat_plot(nn_h_inf_data,"ntransmissions_hyperedge","pairwise")

using CairoMakie
using LaTeXStrings

function alpha_hyperstat_plot(nn_h_inf_data, hyper_stat, hyper_beta_func)
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1],
        xlabel = L"\alpha",
        ylabel = "Hyperedge Size",
        title = "$hyper_stat\n$hyper_beta_func normalization",
        yscale = log10)

    alpha_range = range(0, 2, 15)
    transparency_range = 4
    colors = cgrad(:viridis, length(beta_vals), categorical=true)

    for (ind, beta) in enumerate(beta_vals)
        a_data = reduce(hcat, [nn_h_inf_data[(a, hyper_stat, hyper_beta_func)][:, ind] for a in alpha_range])
        
        lines!(ax, alpha_range, vec(mapslices(x -> nanfindmax(x), a_data, dims=1)),
            # color = (:blue, exp(ind/transparency_range)/exp(length(beta_vals)/transparency_range)),
            color = colors[ind],
            linewidth = 1.8)
    end

    return fig
end

alpha_hyperstat_plot(nn_h_inf_data,"ntransmissions_hyperedge","sqrt")

# alpha_hyperstat_plot(nn_h_inf_data,"ninfected_hyperedge","linear")
# alpha_hyperstat_plot(nn_h_inf_data,"ninfected_hyperedge","sqrt")
# alpha_hyperstat_plot(nn_h_inf_data,"ninfected_hyperedge","pairwise")


# fix beta and look at total infections across alpha 
function total_infections(agg_inf_data,beta,hyper_beta_func)
    ind_to_beta = Dict(enumerate(beta_vals))
    beta_to_ind = Dict(k=>v for (v,k) in pairs(ind_to_beta))

    ind = beta_to_ind[beta]

    a_data = [agg_inf_data[(a,"infections",hyper_beta_func)][ind] for a in range(0,2,15)]
    f = Plots.plot(range(0,2,15),a_data,leg=false,linewidth=2.0,
            xlabel=L"\alpha",
            ylabel="Total Infections",
            title="$hyper_beta_func normalization")
    
    return f 
end

total_infections(agg_inf_data,1e-2,"sqrt")

