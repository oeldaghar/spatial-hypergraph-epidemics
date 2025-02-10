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
            ["spatial-hypergraph-5000-2-0.0-1.jsonl",
                "spatial-hypergraph-5000-2-1.0-1.jsonl",
                "spatial-hypergraph-5000-2-2.0-1.jsonl"],
            # 5k nodes in 5d
            ["spatial-hypergraph-5000-5-0.0-1.jsonl",
                "spatial-hypergraph-5000-5-1.0-1.jsonl",
                "spatial-hypergraph-5000-5-2.0-1.jsonl"],
            # 50k nodes in 2d
            ["spatial-hypergraph-50000-2-0.0-1.jsonl",
                "spatial-hypergraph-50000-2-1.0-1.jsonl",
                "spatial-hypergraph-50000-2-2.0-1.jsonl"],
            # 50k nodes in 5d
            ["spatial-hypergraph-50000-5-0.0-1.jsonl",
                "spatial-hypergraph-50000-5-1.0-1.jsonl",
                "spatial-hypergraph-50000-5-2.0-1.jsonl"],
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

function load_epidemic_data(data_dir,fname;excluded=[])
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
            obj = JSON3.read(line)
            for k in keys(obj)
                key = string(k)
                if !(key in excluded)
                    push!(data[key],obj[key])
                end
            end
        end
    end
    # post processing
    data["infections"] = reduce(hcat,data["infections"])'
    return data 
end

# test_data = load_epidemic_data("data/hysteresis/sirs/", "spatial-hypergraph-5000-2-2.0-1-newalpha-TEST.jsonl")
# data_dir = "data/hysteresis/sirs/"
# fname = "spatial-hypergraph-5000-2-2.0-1-newalpha-TEST.jsonl"
# f1,f2,f3 = plot_hysteresis(data_dir, fname)
# display(f1)
# display(f2)
# display(f3)
# similar to above but rather than load the entire data in, we reduce on the memory intensive keys 
# "ntransmissions_hyperedge" and "ninfected_hyperedge".
function load_epidemic_data_reduce(data_dir,fname,reduce_func;excluded=[])
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
            obj = JSON3.read(line)
            for k in keys(obj)
                key = string(k)
                if !(key in excluded || key in reduce_keys) 
                    push!(data[key],obj[key])
                elseif key in reduce_keys
                    # load it in, reduce and push reduced output to 
                    push!(data[key],reduce_func(obj[key]))
                end
            end
        end
    end
    # post processing
    data["infections"] = reduce(hcat,data["infections"])'
    return data 
end

# data_dir = "data/hysteresis/sirs/"
# fname = "spatial-hypergraph-5000-5-2.0-1_results.jsonl"
# fpath = joinpath(data_dir,fname)
# ispath(fpath)
# data = load_epidemic_data(data_dir,fname,excluded=["ntransmissions_hyperedge","ninfected_hyperedge"])
# data["ntransmissions_hyperedge"] #should return empty 

# plot 1: can we see hysteresis as we vary beta for a single graph 
# plot binned rolling average by starting condition vs beta as scatter plot
# if there is no bistability (not in that parameter regime), then
# we'd expect to see that points corresponding to the same beta should be very close to one another
# (the initial condition doesn't matter)
# otherwise, we'd expect to see meaningfully different values for the same beta 

# Example usage
function plot_hysteresis(data_dir, fname)
    data = load_epidemic_data(data_dir,fname,excluded = ["ntransmissions_hyperedge","ninfected_hyperedge"])

    # group data by starting condition, beta, and hyper_beta_func
    filtered_data = hcat(data["beta"], data["hyper_beta_func"], data["infections"])
    grouped_inds = group_by_function(x -> (x[1], x[2], x[3]), eachrow(filtered_data))
    grouped_rolling_averages = Dict(k => rolling_average(data["infections"][grouped_inds[k], :], 1000) for k in keys(grouped_inds))

    # put the plot together 
    static_keys = sort(collect(keys(grouped_rolling_averages)))
    linear_keys = [k for k in static_keys if k[2] == "linear"]
    sqrt_keys = [k for k in static_keys if k[2] == "sqrt"]

    beta_vals = unique(first.(linear_keys))
    beta_map = Dict(map(x -> (x[2], x[1]), enumerate(sort(beta_vals))))

    # beta_map[beta_val] = x_ind
    function _plot_keys(keys, title_suffix)
        # xs = [beta_map[x[1]] for x in keys]
        xs = [x[1] for x in keys]
        ys = [grouped_rolling_averages[x] for x in keys]

        f = Plots.scatter(xs, ys, leg = false,
                xticks = ([1e-3,1e-2,1e-1,1e0], ["1e-3","1e-2","1e-1","1e0"]),
                xlabel = "Pairwise Infection Probabilihcatty", 
                ylabel = "Binned Average\nRolling Infections", 
                title = "Beta vs Rolling Infections - $title_suffix\n$fname",
                markerstrokewidth = 0,
                xlims=(8e-4,1.05e0))
        Plots.plot!(f,xscale=:log10)
    end

    return _plot_keys(linear_keys, "Linear"),_plot_keys(sqrt_keys, "Sqrt"), _plot_keys([k for k in static_keys if k[2] == "pairwise"], "Pairwise")
end

################################
## Figure 1 - hysteresis testing
################################
# Example usage
# data_dir = "data/hysteresis/sirs/"
# fname = "spatial-hypergraph-50000-5-0.0-1.jsonl"
# f1,f2 = plot_hysteresis(data_dir, fname)
# f1
# f2

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
    mydict = Dict{Any,Float64}()
    for dict in vec_of_dicts
        for key in keys(dict)
            mydict[key] = get(mydict,key,0)+dict[key]
        end
    end
    nitems = lastindex(vec_of_dicts)
    for key in keys(mydict)
        mydict[key] = mydict[key]/nitems
    end
    return mydict
end

function dictionary_rolling_average(vec_of_dicts::S, window::Int)::Dict where S
    # build a dict to track sum prior to average 
    rolling_avgs = Dict{Any, Float64}()
    for dict in vec_of_dicts[end-window:end]
        for (key,val) in pairs(dict)
            rolling_avgs[key] = get(rolling_avgs,key,0)+val
        end
    end
    
    nitems = lastindex(vec_of_dicts)-window+1
    for key in keys(rolling_avgs)
        rolling_avgs[key] = rolling_avgs[key]/nitems
    end
    return rolling_avgs
end

function rolling_hyperedge_stats(data,hyperedge_key="ntransmissions_hyperedge",stream=false)
    # group data by starting condition, beta, and hyper_beta_func
    filtered_data = hcat(data["beta"], data["hyper_beta_func"])
    grouped_inds = group_by_function(x -> (x[1],x[2]), eachrow(filtered_data))
    # compute rolling averages for each trajectory/simulation and average over group size 
    println("Computing rolling averages...")
    newdata = Dict()
    for (key,inds) in pairs(grouped_inds)
        # for original data 
        if !stream
            rolling_avgs = map(ind->dictionary_rolling_average(data[hyperedge_key][ind],1000),inds)
            group_rolling_avg = average_over_dicts(rolling_avgs)
        else
            group_rolling_avg = average_over_dicts(data[hyperedge_key][inds])
        end
        newdata[key] = group_rolling_avg
    end

    # key set up 
    static_keys = collect(keys(newdata))
    linear_keys = [key for key in static_keys if key[2]=="linear"]
    sqrt_keys = [key for key in static_keys if key[2]=="sqrt"]
    hyperedge_keys = string.(keys(newdata[first(keys(newdata))]))
    ymin,ymax = extrema(parse.(Int,hyperedge_keys))
    
    beta_vals = sort(unique(data["beta"]))
    # fs = []
    pdata = Dict()
    for hyper_beta_func in unique(data["hyper_beta_func"])
        # main plotting setup 
        plotting_data = zeros(Float64,ymax,lastindex(beta_vals))
        for (i,beta) in enumerate(beta_vals)
            dict = newdata[(beta,hyper_beta_func)] # filter to linear 
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
        data = load_epidemic_data_reduce(data_dir,fname,x->dictionary_rolling_average(x,1000))
    else
        data = load_epidemic_data(data_dir,fname)    
    end
    # generate individual figures 
    pdata = Dict()
    for hyperkey in ["ntransmissions_hyperedge","ninfected_hyperedge"]
        output = rolling_hyperedge_stats(data,hyperkey,stream)
        pdata[hyperkey] = output 
        # push!(figs,f1)
        # push!(figs,f2)
    end
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
    # "spatial-hypergraph-5000-2",
    # "spatial-hypergraph-5000-5",
    # larger 
    # "spatial-hypergraph-50000-2",
    "spatial-hypergraph-50000-5"
]
for base_graph in base_graph_names
# base_graph = "spatial-hypergraph-50000-2"
    graph_names = map(x->"$(splitext(basename(x))[1]).jsonl",get_fnames(base_graph))

    start_time = time()
    pdata = Dict()
    for fname in graph_names #smaller graph 
        alpha_val = parse(Float64,split(fname,"-")[5])
        
        out = hypergraph_stats_data(data_dir,fname,true)
        pdata[alpha_val] = out
    end
    end_time = time()
    println("PLOTTING DATA GENERATED IN $(end_time-start_time) SECONDS")

    # flatten data
    newdata = Dict()
    for (key1, val1) in pdata
        for (key2, val2) in val1
            for (key3, val3) in val2
                newdata[(key1, key2, key3)] = deepcopy(val3)
            end
        end
    end

    max_size = maximum(size.(values(newdata),1))
    for (key,val) in pairs(newdata)
        newdata[key] = _pad_matrix(val,max_size)
    end


    hypergraph_info = Dict()
    for gname in graph_names #smaller graph 
        fname = get_fnames(splitext(gname)[1])[1]
        
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

    for col_norm in [true, false]
        figs=make_transmission_heatmaps(newdata,hypergraph_info,col_norm)

        # alpha_vals = sort(collect(keys(hypergraph_info)))
        # for alph in alpha_vals
        #     f = figs[(alph, "ntransmissions_hyperedge", "linear")]
        #     Plots.display(f)
        # end
        for key in collect(keys(figs))
            alph,hyperedge_key,normalization = key
            f = figs[key]
            Plots.plot!(f,dpi=1000)

            fig_name = "$hyperedge_key-$base_graph-$alph-$normalization"
            if col_norm
                fig_name = "$fig_name-max_normalized"
            end
            Plots.savefig(f,joinpath(figure_path,"$fig_name.pdf"))
            Plots.savefig(f,joinpath(figure_path,"$fig_name.png"))
        end
    end
end

#= 
keys(figs)
function myplot_touchup(figs)
    fs = deepcopy(figs)
    l = @layout [grid(2, 3, widths=[0.31, 0.31, 0.38])]
    transmission_fig = Plots.plot(
        fs[(0.0,"ntransmissions_hyperedge","linear")],
        fs[(1.0,"ntransmissions_hyperedge","linear")],
        fs[(2.0,"ntransmissions_hyperedge","linear")],
        fs[(0.0,"ntransmissions_hyperedge","sqrt")],
        fs[(1.0,"ntransmissions_hyperedge","sqrt")],
        fs[(2.0,"ntransmissions_hyperedge","sqrt")],
        layout = l,
        size=(2000,800)
    )

    #remove cbar 
    for i in [1,2,4,5]
        Plots.plot!(transmission_fig[i],cbar=false)
    end
    # remove xlabel from first row 
    for i in [1,2,3]
        Plots.plot!(transmission_fig[i],xlabel="")
    end
    # remove ylabels 
    for i in [2,3,5,6]
        Plots.plot!(transmission_fig[i],ylabel="")
    end
    # adjust margins 
    Plots.plot!(transmission_fig[4],bottom_margin=8mm,left_margin=14mm)
    Plots.plot!(transmission_fig[1],top_margin=5mm)
    Plots.plot!(transmission_fig, plot_title="Node-Coverage Normalized Transmissions", plot_titlefontsize=20,
            dpi=2000)
    return transmission_fig
end


myplot_touchup(figs)

figs=make_transmission_heatmaps(newdata,hypergraph_info,true)
myplot_touchup(figs)

figs[(2.0, "ntransmissions_hyperedge", "linear")]
figs[(1.0, "ntransmissions_hyperedge", "linear")]
figs[(0.0, "ntransmissions_hyperedge", "linear")]

Plots.plot(figs[(2.0, "ntransmissions_hyperedge", "linear")],clims=(-3,4))
Plots.plot(figs[(1.0, "ntransmissions_hyperedge", "linear")],clims=(-3,4))
Plots.plot(figs[(0.0, "ntransmissions_hyperedge", "linear")],clims=(-3,4))




##########################
## Individual Trajectories
########################## 
function single_trajectory_hyperedge_transmissions(data,beta=5e-2,hyper_beta_func="linear",initial_num_infected_nodes=2500)
    filtered_data = hcat(data["beta"], data["hyper_beta_func"], data["infections"])
    grouped_inds = group_by_function(x -> (x[1],x[2],x[3]), eachrow(filtered_data))
    # grab the group of inds 
    inds = grouped_inds[(beta,hyper_beta_func,initial_num_infected_nodes)]

    # vec_vec_of_dicts = data["ntransmissions_hyperedge"][inds]

    all_keys = Set()
    for ind in inds 
        for dict in data["ntransmissions_hyperedge"][ind]
            union!(all_keys,keys(dict))
        end
    end
    ymin,ymax = extrema(parse.(Int,all_keys))

    fs = []
    for (i,ind) in enumerate(inds)
        plotting_data = zeros(Float64,ymax,lastindex(data["ntransmissions_hyperedge"][ind]))
        for (timestep,dict) in enumerate(data["ntransmissions_hyperedge"][ind])
            for key in keys(dict)
                plotting_data[parse(Int,key),timestep] = dict[key]
            end
        end
        f = heatmap(plotting_data,clims=(1,20),
                yscale=:log10,
                title = "beta: $beta, norm: $hyper_beta_func, initial infs: $initial_num_infected_nodes, trial: $i",
                xlabel="Time Step",
                ylabel="Transmissions by Hypergraph Size")
        push!(fs,f)
    end
         
    return fs
end

# trajectory figures 
fs = single_trajectory_hyperedge_transmissions(data,1e-1,"sqrt",2500)
for f in fs
    plot!(f,title="$fname\n$(f[1][:title])",
        top_margin=2mm)
end

for i=1:min(lastindex(fs),5)
    # save figure 
    plot!(fs[i],dpi=1000)
    savefig(fs[i],joinpath(figure_path,"hyperedge_transmissions-trajectory_$i-$(splitext(fname)[1]).png"))
    savefig(fs[i],joinpath(figure_path,"hyperedge_transmissions-trajectory_$i-$(splitext(fname)[1]).pdf"))
end

# f = rolling_hyperedge_transmissions(data)
# title!(f,fname)
# plot!(f,dpi=1000)
# savefig(f,joinpath(figure_path,"rolling_hyperedge_transmissions-$(splitext(fname)[1]).png"))
# savefig(f,joinpath(figure_path,"rolling_hyperedge_transmissions-$(splitext(fname)[1]).pdf"))

# fname = "spatial-hypergraph-5000-5-0.0-1_results.jsonl"
# data_dir = "data/hysteresis/sirs/"
# data = load_epidemic_data(data_dir,fname)

# fname
# data


# test on small data set 
data_dir = "data/hysteresis/sirs/"
fname = "spatial-hypergraph-50000-5-1.8571428571428572-1.jsonl"
tmp = load_epidemic_data_reduce(data_dir,fname,x->dictionary_rolling_average(x,1000));
=#