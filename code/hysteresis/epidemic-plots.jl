using StatsPlots
using Plots
using CairoMakie
using Measures
using StatsBase
using JSON 
using JSON3
using ProgressMeter
import LinearAlgebra: norm

println("LOADING IN TOOLS...")
include("../data-io.jl")
include("../hypergraph-tools.jl")

data_dir = "data/hysteresis/sirs/"
figure_path = "data/hysteresis/figures"
if !ispath(figure_path)
    mkpath(figure_path)
end

## GLOBALS 
alpha_vals = collect(range(0,2,15))
beta_vals = vcat(1e-3:1e-3:1e-2, 2e-2:1e-2:1e-1, 2e-1:1e-1:9e-1)
beta_vals = sort(unique(beta_vals))

ind_to_beta = Dict(enumerate(beta_vals))
beta_to_ind = Dict(v=>k for (k,v) in pairs(ind_to_beta))

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


# takes about 3 mins 
aggregated_data = load_aggregated_data("aggregated_data","aggregated-sirs-output.json")
aggregated_data1 = load_aggregated_data("aggregated_data","aggregated-sirs-output-scratch-v1.json")
aggregated_data2 = load_aggregated_data("aggregated_data","aggregated-sirs-output-scratch-v2.json")
aggregated_data3 = load_aggregated_data("aggregated_data","aggregated-sirs-output-scratch-v3.json")

# reduce to graph we care about 
common_keys = reduce(intersect,[keys(aggregated_data),keys(aggregated_data1),keys(aggregated_data2),keys(aggregated_data3)])
common_keys = collect(common_keys)
# look at beta values common to all sims for -50000-2-
b0 = unique(aggregated_data[common_keys[1]]["beta"])
b1 = unique(aggregated_data1[common_keys[1]]["beta"])
b2 = unique(aggregated_data2[common_keys[1]]["beta"])
b3 = unique(aggregated_data3[common_keys[1]]["beta"])
bvals = reduce(intersect,[b0,b1,b2,b3])


d0 = unique(aggregated_data[common_keys[1]]["delta"])
d1 = unique(aggregated_data1[common_keys[1]]["delta"])
d2 = unique(aggregated_data2[common_keys[1]]["delta"])
d3 = unique(aggregated_data3[common_keys[1]]["delta"])
dvals = reduce(union,[d0,d1,d2,d3])


g0 = unique(aggregated_data[common_keys[1]]["gamma"])
g1 = unique(aggregated_data1[common_keys[1]]["gamma"])
g2 = unique(aggregated_data2[common_keys[1]]["gamma"])
g3 = unique(aggregated_data3[common_keys[1]]["gamma"])


# plot 1: can we see hysteresis as we vary beta for a single graph 
# plot binned rolling average by starting condition vs beta as scatter plot
# if there is no bistability (not in that parameter regime), then
# we'd expect to see that points corresponding to the same beta should be very close to one another
# (the initial condition doesn't matter)
# otherwise, we'd expect to see meaningfully different values for the same beta 
function _custom_heatmap(pdata,xvals = -1, max_normalize=false)
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
    if max_normalize
        new_mat = mapslices(x -> x ./ maximum(x[(!).(isnan.(x))]), new_mat, dims=1)
    end

    if xvals == -1
        xvals = 1:lastindex(new_mat,2)
    end
    f = Plots.heatmap(xvals, ybins[1:end-1], log10.(new_mat),
                yscale=:log10,
                color=:viridis,
                )
    return f
end

#### PLOT 1 TOTAL INFECTIONS #####
function total_infections_heatmap(agg_data,gname_key="-50000-2-")
    graph_names = collect(keys(aggregated_data))
    graph_names = filter(x->occursin(gname_key,x),graph_names)
    # sort by value of alpha 
    a_vals = map(x->parse(Float64,split(x,"-")[5]),graph_names)
    p = sortperm(a_vals)
    graph_names = graph_names[p]

    b_vals = unique(agg_data[graph_names[1]]["beta"])
    
    # build data for plot
    linear_data = zeros(Float64,length(b_vals),length(alpha_vals))
    pairwise_data = zeros(Float64,length(b_vals),length(alpha_vals))
    sqrt_data = zeros(Float64,length(b_vals),length(alpha_vals))
    
    for (col,gname) in enumerate(graph_names) # loop over alpha
        data = agg_data[gname]
        # group data 
        filtered_data = hcat(data["beta"], data["hyper_beta_func"])
        grouped_inds = group_by_function(x -> (x[1], x[2]), eachrow(filtered_data))
        grouped_averages = Dict(k => mean(data["infections"][grouped_inds[k],:]) for k in keys(grouped_inds))
        for (row,beta) in enumerate(b_vals)
            linear_data[row,col] = grouped_averages[(beta,"linear")]
            sqrt_data[row,col] = grouped_averages[(beta,"sqrt")]
            pairwise_data[row,col] = grouped_averages[(beta,"pairwise")]
        end
    end    
    return linear_data, sqrt_data, pairwise_data
end

linear_data, sqrt_data, pairwise_data = total_infections_heatmap(aggregated_data,"-50000-2-")


f1 = Plots.heatmap(alpha_vals,beta_vals,log10.(linear_data),yscale=:log10, clims=(0.5,4.4))
Plots.plot!(f1,xlabel = L"\alpha", ylabel = L"\beta",title="Total Infections g(m)=m",
            framestyle=:box,
            thickness_scaling=1.2,
            guidefontsize=15,
            top_margin=0Measures.mm)
Plots.savefig(f1,"data/output/figures/final/total-infections-linear.pdf")

f2 = Plots.heatmap(alpha_vals,beta_vals,log10.(sqrt_data),yscale=:log10, clims=(0.5,4.4))
Plots.plot!(f2,xlabel = L"\alpha", ylabel = L"\beta",title="Total Infections, g(m)=sqrt(m)",
            framestyle=:box,
            thickness_scaling=1.2,
            guidefontsize=15,
            top_margin=0Measures.mm)
Plots.savefig(f2,"data/output/figures/final/total-infections-sqrt.pdf")

f3 = Plots.heatmap(alpha_vals,beta_vals,log10.(pairwise_data),yscale=:log10, clims=(0.5,4.4))
Plots.plot!(f3,xlabel = L"\alpha", ylabel = L"\beta",title="Total Infections, g(m)=1",
            framestyle=:box,
            thickness_scaling=1.2,
            guidefontsize=15,
            top_margin=0Measures.mm)
Plots.savefig(f3,"data/output/figures/final/total-infections-pairwise.pdf")


mytitle = ""
fname_suffix = ""
figs = []
for hyper_beta_func in ["pairwise","sqrt","linear"]
    if hyper_beta_func=="pairwise"
        tmp = pairwise_data
        mytitle = "g(m)=1"
        fname_suffix = "pairwise"
    elseif hyper_beta_func == "sqrt"
        tmp = sqrt_data
        mytitle = "g(m)=sqrt(m)"
        fname_suffix = "sqrt"
    elseif hyper_beta_func == "linear"
        tmp = linear_data
        mytitle = "g(m)=m"
        fname_suffix = "linear"
    end
    f = Plots.plot(leg=false)
    colors = cgrad(:plasma, size(tmp,1), categorical=true)
    for (ind,row) in enumerate(eachrow(tmp))
        Plots.plot!(f,alpha_vals,vec(row),c =colors[ind])
    end
    Plots.plot!(f,ylabel="Total Infections",xlabel=L"\alpha",
                    title=mytitle,
                    title_fontsize=8,
                    framestyle=:box,
                    thickness_scaling=1.2,
                    guidefontsize=15,
                    yscale=:log10)
    display(f)
    Plots.savefig(f,"data/output/figures/final/total-infections-trajectory-$fname_suffix.pdf")
    push!(figs,f)
end

function hysteresis_data(agg_data,gname_key="-50000-2-")
    graph_names = collect(keys(agg_data))
    graph_names = filter(x->occursin(gname_key,x),graph_names)
    # sort by value of alpha 
    a_vals = map(x->parse(Float64,split(x,"-")[5]),graph_names)
    p = sortperm(a_vals)
    graph_names = graph_names[p]
    
    # build data for plot
    linear_data = zeros(Float64,length(beta_vals),length(alpha_vals))
    pairwise_data = zeros(Float64,length(beta_vals),length(alpha_vals))
    sqrt_data = zeros(Float64,length(beta_vals),length(alpha_vals))
    
    function _compute_range(itr)
        min_val,max_val = extrema(itr)
        return max_val - min_val
    end

    for (col,gname) in enumerate(graph_names) # loop over alpha
        data = agg_data[gname]
        # group data 
        filtered_data = hcat(data["beta"], data["hyper_beta_func"])
        grouped_inds = group_by_function(x -> tuple(x...), eachrow(filtered_data))
        for (row,beta) in enumerate(beta_vals)
            linear_data[row,col] = _compute_range(mean(data["infections"][grouped_inds[(beta,"linear")],:],dims=2))
            sqrt_data[row,col] = _compute_range(mean(data["infections"][grouped_inds[(beta,"sqrt")],:],dims=2))
            pairwise_data[row,col] = _compute_range(mean(data["infections"][grouped_inds[(beta,"pairwise")],:],dims=2))
        end
    end    
    return linear_data, sqrt_data, pairwise_data    
end

linear_hys, sqrt_hys, pairwise_hys = hysteresis_data(aggregated_data,"-50000-2-")


mytitle = ""
for hyper_beta_func in ["pairwise", "sqrt", "linear"]
    if hyper_beta_func == "pairwise"
        tmp = pairwise_hys
        mytitle = "Total Infections Range\ng(m)=1"
    elseif hyper_beta_func == "sqrt"
        tmp = sqrt_hys
        mytitle = "Total Infections Range\ng(m)=sqrt(m)"
    elseif hyper_beta_func == "linear"
        tmp = linear_hys
        mytitle = "Total Infections Range\ng(m)=m"
    end
    f = Plots.heatmap(alpha_vals, beta_vals, tmp,
            yscale=:log10,
            clims=(0,600),
            xlabel = L"\alpha",
            ylabel = L"\beta",
            guidefontsize=15,
            title = mytitle,
            top_margin=3Plots.mm,
            left_margin=-2Plots.mm,
            bottom_margin=-2Plots.mm,
            thickness_scaling = 1.2,
            size=(600,450)
        )
    Plots.savefig(f,"data/output/figures/final/hysteresis-$hyper_beta_func.pdf")
    display(f)
end



# Example usage
function plot_hysteresis(agg_data::Dict,fname::String)    
    data = agg_data[fname]
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
        # data specific
        ymax = occursin("-50000-",fname) ? 30000 : 3000
        if maximum(ys) >= ymax
            @warn("ymax smaller than observed max, data may have changed")
        end
        parts = split(fname,"-")
        truncated_title = "n: $(parts[3]), d: $(parts[4]), alpha: $(round(parse(Float64,parts[5]),digits=4))"
        f = Plots.scatter(xs, ys, leg = false,
                xticks = ([1e-3,1e-2,1e-1,1e0], ["1e-3","1e-2","1e-1","1e0"]),
                xlabel = "Pairwise Infection Probability", 
                ylabel = "Binned Average\nRolling Infections", 
                title = "Beta vs Rolling Infections - $title_suffix\n$truncated_title",
                markerstrokewidth = 0,
                xlims=(8e-4,1.05e0),
                ylims=(-750,ymax)
                )
        Plots.plot!(f,xscale=:log10)
    end

    return _plot_keys(linear_keys, "Linear"),_plot_keys(sqrt_keys, "Sqrt"), _plot_keys(pairwise_keys, "Pairwise")
end

################################
## Figure 1 - hysteresis testing
################################



# fname = "spatial-hypergraph-50000-5-0.8571428571428571-1-newalpha.jsonl"

figure_path = "data/hysteresis/figures/binned-infections/"
if !ispath(figure_path)
    mkpath(figure_path)
end
@showprogress for key in keys(aggregated_data)
    f1,f2,f3 = plot_hysteresis(aggregated_data,key)
    # save figures 
    Plots.plot!(f1,dpi=300)
    Plots.plot!(f2,dpi=300)
    Plots.plot!(f3,dpi=300)
    gname = splitext(key)[1]
    Plots.savefig(f1, joinpath(figure_path,"linear-$gname.pdf"))
    Plots.savefig(f1, joinpath(figure_path,"linear-$gname.png"))

    Plots.savefig(f2, joinpath(figure_path,"sqrt-$gname.pdf"))
    Plots.savefig(f2, joinpath(figure_path,"sqrt-$gname.png"))

    Plots.savefig(f3, joinpath(figure_path,"pairwise-$gname.pdf"))
    Plots.savefig(f3, joinpath(figure_path,"pairwise-$gname.png"))
end

        
#############################
## PLOTS for hypergraph stats
#############################

## TESTING - doesn't make it into the paper
for fname in graph_names
    data = aggregated_data[fname]

    filtered_data = hcat(data["beta"], data["hyper_beta_func"],data["initial_num_infected_nodes"])
    grouped_inds = group_by_function(x -> (x[1], x[2]), eachrow(filtered_data))

    max_norm = 0
    item = nothing
    for (key,group) in pairs(grouped_inds)
        tmp = data["ntransmissions_hyperedge"][group,:]'
        qinfo = mapslices(x->quantile(x,[0,0.5,1]),tmp,dims=2)
        rel_norm = norm(qinfo[:,3].-qinfo[:,1])/norm(qinfo[:,1])
        if rel_norm>max_norm
            max_norm = rel_norm
            item = key
        end
    end
    println("group: $key relative norm: $(rel_norm)")
end

function plot_with_ribbons(qinfo)
    x = 1:size(qinfo, 1)
    lower = qinfo[:, 1]
    middle = qinfo[:, 2]
    upper = qinfo[:, 3]

    Plots.plot(x, middle, 
         ribbon = (middle - lower, upper - middle),
         fillalpha = 0.3,
         linewidth = 2,
         color = :blue,
         xlabel = "Hyperedge Size",
         ylabel = "Average Transmissions",
         title = "",
         legend = false,
         size = (800, 600))
end

# maximum norm 
# (0.002, "sqrt") relative norm: 8.894601532970313)
# inds = grouped_inds[(2e-3,"sqrt")]
# tmp = data["ntransmissions_hyperedge"][inds,:]'
# qinfo = mapslices(x->quantile(x,[0,0.5,1]),tmp,dims=2)
# fig = plot_with_ribbons(qinfo)
# # can ignore, these are below epidemic thresholds
### CONCLUSION 1
#  AVERAGED HYPEREDGE TRANSMISSIONS STABLE ACROSS TRAJECTORIES 
#  BETA + HYPER_BETA_FUNC enough


# bin the data into coarser bins 
# ind_to_beta[5]
# Plots.plot(sqrt_data[5,:],yscale=:log10)
# graph_names = collect(keys(aggregated_data))
# hinfo_dict = Dict()
# graph_data_dir = "data/hypergraphs/"
# @showprogress for graph_jsonl in graph_names
#     hinfo_dict[graph_jsonl] = Dict()
#     gname = "$(first(splitext(graph_jsonl))).txt"
#     # load graph 
#     node_coverage,edge_coverage = hyperedge_information(joinpath(graph_data_dir,gname))
#     hinfo_dict[graph_jsonl]["node_coverage"] = node_coverage
#     hinfo_dict[graph_jsonl]["edge_coverage"] = edge_coverage
# end

function get_node_coverage(aggregated_data::Dict,hinfo_dict::Dict,graph_name="-50000-2-")
    graph_names = collect(keys(aggregated_data))
    graph_names = filter(x->occursin(graph_name,x),graph_names)
    # sort by value of alpha 
    a_vals = map(x->parse(Float64,split(x,"-")[5]),graph_names)
    p = sortperm(a_vals)
    graph_names = graph_names[p]

    max_hsize = maximum(map(x->size(aggregated_data[x]["ntransmissions_hyperedge"],2),graph_names))
    node_coverage = zeros(Float64,max_hsize,length(alpha_vals))
    for (col,gname) in enumerate(graph_names)
        tmp = hinfo_dict[gname]["node_coverage"]
        node_coverage[:,col] = vcat(tmp,zeros(Float64,max_hsize-length(tmp)))
    end
    return node_coverage
end

function get_beta_transmissions_data(agg_data,graph_name="-50000-2-",beta=2e-2,hyper_beta_func="linear";avg_inf_th=500,coverage=false)
    # get related keys 
    graph_names = collect(keys(agg_data))
    graph_names = filter(x->occursin(graph_name,x),graph_names)
    # sort by value of alpha 
    a_vals = map(x->parse(Float64,split(x,"-")[5]),graph_names)
    p = sortperm(a_vals)
    graph_names = graph_names[p]

    # get maximum h_size over all keys and initialize plotting data 
    max_hsize = maximum(map(x->size(agg_data[x]["ntransmissions_hyperedge"],2),graph_names))
    pdata = zeros(Float64,max_hsize,length(alpha_vals))
    # populate data 
    for (col,gname) in enumerate(graph_names)
        data = agg_data[gname]
        # group epidemics and pick out our group
        filtered_data = hcat(data["beta"], data["hyper_beta_func"])
        grouped_inds = group_by_function(x -> (x[1], x[2]), eachrow(filtered_data))
        inds = grouped_inds[(beta,hyper_beta_func)]
        # look at infections 
        avg_infs = mean(data["infections"][inds,:])
        if avg_infs>avg_inf_th
        # populate data 
            target = vec(mean(data["ntransmissions_hyperedge"][inds,:],dims=1))
            for (row,val) in enumerate(target)
                pdata[row,col] = val
            end
        end
        # normalize by coverage 
        if coverage
            node_coverage = hinfo_dict[gname]["node_coverage"]
            node_coverage = vcat(node_coverage,zeros(Float64,size(pdata,1)-length(node_coverage)))
            pdata[:,col] ./= node_coverage 
        end
    end
    return pdata
end

# transmissions without node coverage
hyper_beta = "pairwise"
beta = 2e-3
pdata = get_beta_transmissions_data(aggregated_data,
                                "-50000-2-",
                                # 2e-3,
                                beta,
                                hyper_beta,
                                avg_inf_th=0,
                                coverage = false,
                            )
pdata = mapslices(x->x./sum(x), pdata,dims=1)
f = _custom_heatmap(pdata,alpha_vals,false)
Plots.plot!(f,
                xlabel=L"\alpha",
                ylabel="Hyperedge Size",
                title="Normalized Transmissions\nβ=$beta, g(m)=1",
                # clims=(-2.8,0),
                right_margin=3Plots.mm,
                top_margin=3Plots.mm,
                bottom_margin=-1Plots.mm,
                thickness_scaling=1.2)

Plots.savefig(f,"data/output/figures/final/transmissions-2e-3.pdf")
    
# total infections 
f = Plots.scatter(alpha_vals, pairwise_data[beta_to_ind[2e-3],:],
                    markerstrokewidth = 0,
                    markersize=6,
                    yscale=:log10, 
                    leg=false,
                    xlabel = L"\alpha",
                    ylabel="Total Infections",
                    title = L"\beta=0.002, g(m) = 1", 
                    ylims = (1,25000),
                    thickness_scaling = 1.3,
                )
Plots.savefig(f,"data/output/figures/final/total-infections-50000-2-beta_2e-3-pairwise.pdf")


################## TRAILING INFECTIONS ######################

# TODO finish writing this function and compute trailing infections with error bands 
function get_inf_data(agg_data,gname_key="-50000-2-",beta=1e-1,hyper_beta_func="linear")
    # get related keys 
    graph_names = collect(keys(agg_data))
    graph_names = filter(x->occursin(gname_key,x),graph_names)
    # sort by value of alpha 
    a_vals = map(x->parse(Float64,split(x,"-")[5]),graph_names)
    p = sortperm(a_vals)
    graph_names = graph_names[p]

    result = Vector{Vector}()
    for gname in graph_names 
        data = agg_data[gname]
        filtered_data = hcat(data["beta"], data["hyper_beta_func"])
        grouped_inds = group_by_function(x->(x[1],x[2]),eachrow(filtered_data))

        inds = grouped_inds[(beta,hyper_beta_func)]
        push!(result, map(x->mean(data["infections"][x,:]), inds))
    end
    return reduce(hcat,result)
end

# should be same size but may not be 
figs = []
for hyper_beta_func in ["linear","pairwise","sqrt"]
    for b in [1e-1,4e-1,9e-1]
        tmp = get_inf_data(aggregated_data,"-50000-2-",b,hyper_beta_func)
        qinfo = mapslices(x->quantile(x,[0.0,0.5,1.0]),tmp,dims=1)'
        f = Plots.plot(alpha_vals,qinfo[:,2],
                ribbon = (qinfo[:,2] - qinfo[:,1], qinfo[:,3] - qinfo[:,2]),
                fillalpha = 0.3,
                linewidth = 2,
                linestyle = :dash,
                color = :blue,
                xlabel = L"\alpha",
                ylabel = "Total Infections",
                title = "n=50000, d=2\nbeta=$b, g(m)=$hyper_beta_func",
                legend = false,
        )
        Plots.scatter!(f,alpha_vals,qinfo[:,2],color = :blue,
                    markerstrokewidth = 0,
                    markersize=5,
        )
        display(f)
        push!(figs,f)
    end
end
Plots.plot(figs...,layout=(3,3),size=(900,900))