# initial figures for hysteresis effects 
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
# aggregated_data = load_aggregated_data("aggregated-sirs-output.json")
# aggregated_data = load_aggregated_data("aggregated-sirs-output-scratch-v1.json")
aggregated_data = load_aggregated_data("aggregated-sirs-output.json")

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
    if max_normalize
        new_mat = mapslices(x -> x ./ maximum(x[(!).(isnan.(x))]), new_mat, dims=1)
    end
    # yvals = [mean(ybins[i:i+1]) for i=1:lastindex(ybins)-1]
    # xrange, yrange, data_matrix, yaxis scale
    if xvals == -1
        xvals = 1:lastindex(new_mat,2)
    end
    f = Plots.heatmap(xvals, ybins[1:end-1], log10.(new_mat),
    # f = Plots.heatmap(1:lastindex(new_mat,2), yvals, log10.(new_mat),
                yscale=:log10,
                color=:viridis,
                # clims=(-5,0),
                )
    return f#,new_mat,ybins 
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

Plots.plot(figs...,layout=(1,3),size=(1200,300),
        top_margin=3Plots.mm,
        bottom_margin = 5Plots.mm,
        left_margin=3Plots.mm,
        link=:y,)





# touch things up if needed and save them

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

for fname in graph_names
    fname = "spatial-hypergraph-50000-2-0.8571428571428571-1-newalpha.jsonl"
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

inds = grouped_inds[(2e-3,"sqrt")]
tmp = data["ntransmissions_hyperedge"][inds,:]'
qinfo = mapslices(x->quantile(x,[0,0.5,1]),tmp,dims=2)

fig = plot_with_ribbons(qinfo)
# can ignore, these are below epidemic thresholds

### CONCLUSION 1
#  AVERAGED HYPEREDGE TRANSMISSIONS STABLE ACROSS TRAJECTORIES 
#  BETA + HYPER_BETA_FUNC enough


# bin the data into coarser bins 
ind_to_beta[5]
Plots.plot(sqrt_data[5,:],yscale=:log10)

graph_names = collect(keys(aggregated_data))

hinfo_dict = Dict()
graph_data_dir = "data/hypergraphs/"
@showprogress for graph_jsonl in graph_names
    hinfo_dict[graph_jsonl] = Dict()
    gname = "$(first(splitext(graph_jsonl))).txt"
    # load graph 
    node_coverage,edge_coverage = hyperedge_information(joinpath(graph_data_dir,gname))
    hinfo_dict[graph_jsonl]["node_coverage"] = node_coverage
    hinfo_dict[graph_jsonl]["edge_coverage"] = edge_coverage
end


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

### FIRST SET OF FIGURES ####
# node_cov = get_node_coverage(aggregated_data, hinfo_dict,"-50000-2-")
# f = _custom_heatmap(node_cov,alpha_vals)
# Plots.plot!(f,xlabel=L"\alpha",ylabel="Hyperedge Size",title="Normalized Node Coverage\nn=50000, d=2",
#             clims=(-5.0,0),
#             top_margin=3Plots.mm)
# Plots.savefig(f,"data/output/figures/final/normalized_node_coverage-50000-2.pdf")

# pdata = get_beta_transmissions_data(aggregated_data,
#                                 "-50000-2-",
#                                 2e-3,
#                                 "pairwise",
#                                 avg_inf_th=0,
#                                 coverage = true,
#                             )
# f = _custom_heatmap(pdata,alpha_vals)
# Plots.plot!(f,
#                 xlabel=L"\alpha",
#                 ylabel="Hyperedge Size",
#                 title="Node Normalized Transmissions\nn=50000, d=2, g(m)=1",
#                 clims=(-5.0,0),
#                 right_margin=3Plots.mm,
#                 top_margin=3Plots.mm)
# Plots.savefig(f,"data/output/figures/final/node_normalized_transmissions-50000-2_beta_2e-3.pdf")


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
    

# tmp = mapslices(c->c./sum(c),pdata,dims=1)
# tmp2 = tmp[end:-1:1,:]
# tmp3 = mapslices(x->cumsum(x),tmp2,dims=1) 
# cols = tmp3[:,[2,5,8,11,14]]
# f = Plots.plot(yscale=:log10,leg=:topleft) # create a new plot
# for col in eachcol(cols)
#     ind = findfirst(col .> 0)
#     Plots.plot!(ind:length(col), col[ind:end])
# end
# f

plot_with_ribbons
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

## END FIRST SET OF FIGURES 

#### SECOND SET OF PLOTS - ROUGHLY SAME TOTAL INFS ####
# total infections 
for b in [2e-2] #5e-1, 6e-1, 7e-1, 8e-1, 9e-1]
    for hyper_beta_func in ["linear", "pairwise", "sqrt"]
        if hyper_beta_func == "pairwise"
            ys = pairwise_data[beta_to_ind[b],:]
            mytitle = "beta=$b, g(m) = 1"
        elseif hyper_beta_func == "linear"
            ys = linear_data[beta_to_ind[b],:]
            mytitle = "beta=$b, g(m) = m"
        elseif hyper_beta_func == "sqrt"
            ys = sqrt_data[beta_to_ind[b],:]
            mytitle = "beta=$b, g(m) = sqrt{m}"
        end
        f = Plots.scatter(alpha_vals, ys,
                        markerstrokewidth = 0,
                        markersize=6,
                        leg=false,
                        xlabel = L"\alpha",
                        ylabel="Total Infections",
                        title = mytitle,
                        ylims = (11000,24000),
                        thickness_scaling=1.3
                    )    
        display(f)
        # Plots.savefig(f,"data/output/figures/final/total-infections-50000-2-beta_4e-1-$hyper_beta_func.pdf")
    end
end




graph_names = collect(keys(aggregated_data))
graph_names = filter(x->occursin("-50000-2-",x),graph_names)
# sort by value of alpha 
a_vals = map(x->parse(Float64,split(x,"-")[5]),graph_names)
p = sortperm(a_vals)
graph_names = graph_names[p]


include("../sirs.jl")

tmp_gname = graph_names[end]
tmp_gname = "$(splitext(tmp_gname)[1]).txt"
hedges = read_hedges("data/hypergraphs/$tmp_gname")
res = first(sirs_hyper(hedges,sample(1:50000,500),5e-1,1e-1,1e-2,5/1000000,1000,"pairwise"))

Plots.plot(res)

println(mean(res[end-100,end]))


# heatmaps 
cov = false
for hyper_beta_func in ["pairwise", "sqrt", "linear"]
    b = 4e-1
    pdata = get_beta_transmissions_data(aggregated_data,
                                "-50000-2-",
                                b,
                                hyper_beta_func,
                                avg_inf_th=0,
                                coverage = cov,
                            )
    f = _custom_heatmap(pdata,alpha_vals)
    if cov 
        mytitle = "Node Normalized Transmissions, β=$b"
    else
        mytitle = "Normalized Transmissions, β=$b"
    end
    if hyper_beta_func=="linear"
        other_stats = "n=50000, d=2, g(m)=m"
    elseif hyper_beta_func == "sqrt"
        other_stats = "n=50000, d=2, g(m)=sqrt(m)"
    elseif hyper_beta_func == "pairwise"
        other_stats = "n=50000, d=2, g(m)=1"
    end
    mytitle = mytitle*"\n"*other_stats
    Plots.plot!(f,
                    xlabel=L"\alpha",
                    ylabel="Hyperedge Size",
                    title=mytitle,
                    clims=(-5.0,0),
                    right_margin=3Plots.mm,
                    top_margin=3Plots.mm)
    display(f)
    # Plots.savefig(f,"data/output/figures/final/transmissions-50000-2_beta_4e-1-$hyper_beta_func.pdf")
end



### BETA PLOTS 
graph_names = collect(keys(aggregated_data))
graph_names = filter(x->occursin("-50000-2-",x),graph_names)
# sort by value of alpha 
a_vals = map(x->parse(Float64,split(x,"-")[5]),graph_names)
p = sortperm(a_vals)
graph_names = graph_names[p]



function get_alpha_transmissions_data(aggregated_data,graph_name="-50000-2-",alpha=0.0,hyper_beta_func="linear";avg_inf_th=0,coverage=false)
    # get related keys 
    graph_names = collect(keys(aggregated_data))
    graph_names = filter(x->occursin(graph_name,x),graph_names)
    # sort by value of alpha 
    a_vals = map(x->parse(Float64,split(x,"-")[5]),graph_names)
    ind = findfirst(a_vals.==alpha)
    
    gname = graph_names[ind]
    
    # get maximum h_size over all keys and initialize plotting data 
    max_hsize = maximum(map(x->size(aggregated_data[x]["ntransmissions_hyperedge"],2),graph_names))
    pdata = zeros(Float64,max_hsize,length(beta_vals))
    # populate data 
    data = aggregated_data[gname]
    # group epidemics and pick out our group
    filtered_data = hcat(data["beta"], data["hyper_beta_func"])
    grouped_inds = group_by_function(x -> (x[1], x[2]), eachrow(filtered_data))
    for (col,beta) in enumerate(beta_vals)
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

tmp = get_alpha_transmissions_data(aggregated_data,"-50000-2-",
                                        alpha_vals[3],
                                        "sqrt",
                                        coverage=true)
f = _custom_heatmap(tmp, beta_vals)
Plots.plot!(f,
                xscale=:log10,
                xlabel=L"\beta",
                ylabel="Hyperedge Size",
                title="Normalized Transmissions")




################## TRAILING INFECTIONS ######################
aggregated_data


for ind = 1:27
    f = Plots.scatter(alpha_vals, linear_data[ind,:],label="g(m)=m",
                        markerstrokewidth=0, markersize = 5,
                    )
    Plots.scatter!(f, alpha_vals, sqrt_data[ind,:],label="g(m)=sqrt(m)",
                        markerstrokewidth=0, markersize = 5,
                    )
    Plots.scatter!(f, alpha_vals, pairwise_data[ind,:],label="g(m)=1",
                        markerstrokewidth=0, markersize = 5,
                    )
    display(f)
end

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




             













# PLOT 1 : hyperedge with most transmissions 
# ignore values where infections not large enough 
# point: alpha + beta can impact which hyperedge edges dominate spreading
# beta large then smaller hyperedges spread faster than larger 
# beta small, then other way around 

ntransmissions_mats = map(beta->get_beta_transmissions_data(aggregated_data,"-50000-2-",beta,"linear",avg_inf_th=0),beta_vals)

for mat in ntransmissions_mats


end


# at threshold, large hyperedges domiante transmissions (superspreading, localization, etc)

# ANOTHER PLOT 
# ventillation biases towards smaller hyperedges but may not move
# total infections depending on parameters 

# ANOTHER PLOT : normalized hyperedge with most transmissions 
# ignore values where infections not large enough 






























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





