# initial figures for hysteresis effects 
using StatsPlots
using Plots
using Measures
using StatsBase
using JSON 

data_dir = "data/hysteresis/sirs/"
figure_path = "data/hysteresis/figures"
if !ispath(figure_path)
    mkpath(figure_path)
end
# gname = "spatial-hypergraph-5000-2"
# fnames = filter(x->occursin(gname,x),readdir(data_dir))

# alphas = map(x->parse(Float64,split(x,"-")[5]),fnames)
# perm = sortperm(alphas)
# alphas = alphas[perm]
# fnames = fnames[perm]


original_graph_names = [
            ["spatial-hypergraph-5000-2-0.0-1-counter_1.json",
                "spatial-hypergraph-5000-2-1.0-1-counter_1.json",
                "spatial-hypergraph-5000-2-2.0-1-counter_1.json"],
            #50k node graphs 
            ["spatial-hypergraph-50000-2-0.0-1.json",
                "spatial-hypergraph-50000-2-1.0-1.json",
                "spatial-hypergraph-50000-2-2.0-1.json"],
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

function load_epidemic_data(data_dir,fname)
    fpath = joinpath(data_dir, fname)
    println("READING DATA...$fpath")
    data = JSON.parsefile(fpath)
    println("EXTRACTING INFO AND CONVERTING TYPES...$fpath")
    data = Dict{String, Array}(data)
    data["ntransmissions_hyperedge"] = Vector{Vector{Dict{String,Int}}}(data["ntransmissions_hyperedge"])
    data["ninfected_hyperedge"] = Vector{Vector{Dict{String,Int}}}(data["ninfected_hyperedge"])
    data["beta"] = Float64.(data["beta"])
    data["hyper_beta_func"] = String.(data["hyper_beta_func"])
    data["infections"] = hcat(Vector{Int}.(data["infections"])...)'
    return data
end
# plot 1: can we see hysteresis as we vary beta for a single graph 
# plot binned rolling average by starting condition vs beta as scatter plot
# if there is no bistability (not in that parameter regime), then
# we'd expect to see that points corresponding to the same beta should be very close to one another
# (the initial condition doesn't matter)
# otherwise, we'd expect to see meaningfully different values for the same beta 

# Example usage
function plot_hysteresis(data_dir, fname)
    # alpha_val = Float64(split(fname,"-")[5])
    data = load_epidemic_data(data_dir,fname)
    # fpath = joinpath(data_dir, fname)
    # println("READING DATA...$fpath")
    # data = JSON.parsefile(fpath)
    # println("EXTRACTING INFO AND CONVERTING TYPES...$fpath")
    # data = Dict{String, Array}(key => data[key] for key in ["beta", "hyper_beta_func", "infections"])
    # data["beta"] = Float64.(data["beta"])
    # data["hyper_beta_func"] = String.(data["hyper_beta_func"])
    # data["infections"] = hcat(Vector{Int}.(data["infections"])...)'

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
    function plot_keys(keys, title_suffix)
        # xs = [beta_map[x[1]] for x in keys]
        xs = [x[1] for x in keys]
        ys = [grouped_rolling_averages[x] for x in keys]

        f = scatter(xs, ys, leg = false,
                xticks = ([1e-3,1e-2,1e-1,1e0], ["1e-3","1e-2","1e-1","1e0"]),
                xlabel = "Pairwise Infection Probabilihcatty", 
                ylabel = "Binned Average\nRolling Infections", 
                title = "Beta vs Rolling Infections - $title_suffix\n$fname",
                markerstrokewidth = 0,
                xlims=(8e-4,1.05e0))
        plot!(f,xscale=:log10)
    end

    return plot_keys(linear_keys, "Linear"),plot_keys(sqrt_keys, "Sqrt")
end

################################
## Figure 1 - hysteresis testing
################################
# Example usage
# data_dir = "data/hysteresis/sirs/"
# fname = "spatial-hypergraph-5000-2-0.0-1.json"
# f1,f2 = plot_hysteresis(data_dir, fname)
# f1

fs = []
for fname in [
            "spatial-hypergraph-5000-2-0.0-1-counter_1.json",
            "spatial-hypergraph-5000-2-1.0-1-counter_1.json",
            "spatial-hypergraph-5000-2-2.0-1-counter_1.json"
        ]
    println("WORKING ON $fname")
    f1,f2 = plot_hysteresis(data_dir, fname)
    ylims!(f1,(-50,3000))
    ylims!(f2,(-50,3000))
    push!(fs,f1)
    push!(fs,f2)
end

# side by side plot as we vary alpha across columns
f_main = plot(fs[1], fs[3], fs[5], fs[2], fs[4], fs[6],
              layout=(2, 3),
              size=(1000, 600),
              left_margin=5mm, right_margin=2mm, top_margin=2mm, bottom_margin=5mm,
              title=["alpha=0" "alpha=1" "alpha=2" "alpha=0" "alpha=1" "alpha=2"],
              titlefontsize=12,
              plot_title="Linear Normalization",
              plot_titlefontsize=18)
# remove xlabel, xtick labels, and adjust spacing for first row 
for i=1:3
    plot!(f_main[i], xticks=(xticks(f_main[i])[1], ""),
        xlabel="",
        bottom_margin=8mm,
        top_margin=8mm)
end
# remove ylabel for 2nd column and beyond
for i=1:6
    if i%3!=1
        plot!(f_main[i], yticks=(yticks(f_main[i])[1], ""),
            ylabel="",
            left_margin=3mm)
    end
end
# Add row titles using annotation
plot!(f_main, annotation=(5e-2, 3800, text("Square Root Normalization", 18, :center)),
                subplot=5)
# save figure 
savefig(f_main, joinpath(figure_path,"hysteresis_plot_v1.pdf"))
plot!(f_main,dpi=1000)
savefig(f_main, joinpath(figure_path,"hysteresis_plot_v1.png"))

#zooming in on one piece 
function stochastic_evolution(data_dir,fname)
    load_epidemic_data(data_dir,fname)
    # fpath = joinpath(data_dir, fname)
    # println("READING DATA...$fpath")
    # data = JSON.parsefile(fpath)
    # println("EXTRACTING INFO AND CONVERTING TYPES...$fpath")
    # data = Dict{String, Array}(key => data[key] for key in ["beta", "hyper_beta_func", "infections"])
    # data["beta"] = Float64.(data["beta"])
    # data["hyper_beta_func"] = String.(data["hyper_beta_func"])
    # data["infections"] = hcat(Vector{Int}.(data["infections"])...)'

    # group data by starting condition, beta, and hyper_beta_func
    filtered_data = hcat(data["beta"], data["hyper_beta_func"], data["infections"])
    grouped_inds = group_by_function(x -> (x[1], x[2], x[3]), eachrow(filtered_data))

    # put the plot together 
    static_keys = sort(collect(keys(grouped_inds)))
    linear_keys = [k for k in static_keys if k[2] == "linear"]
    sqrt_keys = [k for k in static_keys if k[2] == "sqrt"]

    beta_vals = unique(first.(linear_keys))
    beta_map = Dict(map(x -> (x[2], x[1]), enumerate(sort(beta_vals))))

    beta = 1e-2
    # pick a value of beta and do a plot for each normalization
    inds = (first.(static_keys).==beta) .& ([x[2]=="linear" for x in static_keys])

    f = plot()
    for (i,key) in enumerate(static_keys[inds])
        group = grouped_inds[key]
        plot!(f,data["infections"][group,1:200]',c=i,leg=false)
    end

    alpha_val = split(fname,"-")[5]
    plot!(f,
        xlabel="Time Step",
        ylabel="Total Infections",
        title = "alpha=$alpha_val, beta = $beta\$fname")
    return f 
end

filtered_fnames = [
            "spatial-hypergraph-5000-2-0.0-1-counter_1.json",
            "spatial-hypergraph-5000-2-1.0-1-counter_1.json",
            "spatial-hypergraph-5000-2-2.0-1-counter_1.json"
        ]
fs = []
for fname in filtered_fnames
    f = stochastic_evolution(data_dir,fname)
    ylims!(f,(0,5000))
    push!(fs,f)
end

#put them together 
main_f = plot(fs...,layout=(3,1),
        size = (250,650),
        # margins=8mm,
        left_margin = 7mm,
        titlefontsize=11)
savefig(main_f, joinpath(figure_path,"hysteresis_stochastic_plot_v1.pdf"))
plot!(main_f,dpi=1000)
savefig(main_f, joinpath(figure_path,"hysteresis_stochastic_plot_v1.png"))
        
#############################
## PLOTS for hypergraph stats 
#############################
# average over dicts 
function average_over_dicts(vec_of_dicts::Vector{Dict{S,T}}) where {S,T}
    mydict = Dict{S,Float64}()
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

function dictionary_rolling_average(vec_of_dicts::Vector{Dict{S,T}}, window::Int) where {S,T}
    # build a dict to track sum prior to average 
    rolling_avgs = Dict{S,Float64}()
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

function rolling_hyperedge_stats(data,hyperedge_key="ntransmissions_hyperedge")
    # group data by starting condition, beta, and hyper_beta_func
    filtered_data = hcat(data["beta"], data["hyper_beta_func"])
    grouped_inds = group_by_function(x -> (x[1],x[2]), eachrow(filtered_data))
    # compute rolling averages for each trajectory/simulation and average over group size 
    println("Computing rolling averages...")
    newdata = Dict()
    for (key,inds) in pairs(grouped_inds)
        rolling_avgs = map(ind->dictionary_rolling_average(data[hyperedge_key][ind],1000),inds)
        group_rolling_avg = average_over_dicts(rolling_avgs)
        newdata[key] = group_rolling_avg
    end

    # key set up 
    static_keys = collect(keys(newdata))
    linear_keys = [key for key in static_keys if key[2]=="linear"]
    sqrt_keys = [key for key in static_keys if key[2]=="sqrt"]
    hyperedge_keys = keys(newdata[first(keys(newdata))])
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
                plotting_data[parse(Int,key),i] = dict[key]
            end
        end

        pdata[hyper_beta_func] = plotting_data
    end
    return pdata 
end

##################
## Hyperedge stats 
##################
# fname = "spatial-hypergraph-5000-2-0.0-1-counter_1.json"
# fpath = joinpath(data_dir, fname)
# data = load_epidemic_data(data_dir,fname)

function hypergraph_stats_fig(data_dir,fname)
    # load data 
    data = load_epidemic_data(data_dir,fname)    
    # generate individual figures 
    pdata = Dict()
    for hyperkey in ["ntransmissions_hyperedge","ninfected_hyperedge"]
        output = rolling_hyperedge_stats(data,hyperkey)
        pdata[hyperkey] = output 
        # push!(figs,f1)
        # push!(figs,f2)
    end
    return pdata
end


start_time = time()
pdata = Dict()
for fname in original_graph_names[1] #smaller graph 
    alpha_val = parse(Float64,split(fname,"-")[5])
    
    out = hypergraph_stats_fig(data_dir,fname)
    pdata[alpha_val] = out
end
end_time = time()
println("ORIGNAL FIGURES COMPUTED IN $(end_time-start_time) SECONDS")

function _pad_matrix(matrix, nrows=100)
    nnewrows = nrows - size(matrix, 1)
    ncols = size(matrix, 2)
    if nnewrows > 0
        matrix = vcat(matrix, zeros(nnewrows, ncols))
    end
    return matrix
end

# 
newdata = Dict()
for (key1, val1) in pdata
    for (key2, val2) in val1
        for (key3, val3) in val2
            newdata[(key1, key2, key3)] = val3
        end
    end
end

max_size = maximum(size.(values(newdata),1))
for (key,val) in pairs(newdata)
    newdata[key] = _pad_matrix(val,max_size)
end

# plot individual figures 
function make_transmission_heatmaps(data_dict)
    fig_dict = Dict()
    for (key,plotting_data) in pairs(data_dict)
        alpha,hyperkey,hyper_beta_func = key

        beta_vals = vcat(1e-3:1e-3:1e-2, 2e-2:1e-2:1e-1, 2e-1:1e-1:9e-1)
        beta_tick_map = Dict(reverse.(enumerate(beta_vals)))
        ticks = (vcat([beta_tick_map[x] for x in [1e-3,1e-2,1e-1]],lastindex(beta_vals)) ,["1e-3","1e-2","1e-1","1e0"])
        if hyperkey == "ntransmissions_hyperedge"
            mylabel = "Rolling Average Transmission\n by Hyperedge Size"
        elseif hyperkey == "ninfected_hyperedge"
            mylabel = "Rolling Average Infection\n by Hyperedge Size"
        end
        
        f = heatmap(1 .+log10.(plotting_data),
                yscale=:log10,            
                # clims = (-1,20),
                cmap=:viridis,
                ylims=(0.9,size(plotting_data,1)+10),
                xticks = ticks,
                xlabel = "Pairwise Infection Probability",
                ylabel = mylabel,
                title="Alpha:$alpha Norm: $hyper_beta_func")
        fig_dict[key] = f
    end
    return fig_dict
end

figs=make_transmission_heatmaps(newdata)

# put side by side 
l = @layout [grid(2, 3, widths=[0.31, 0.31, 0.38])]
transmission_fig = plot(
    figs[(0.0,"ntransmissions_hyperedge","linear")],
    figs[(1.0,"ntransmissions_hyperedge","linear")],
    figs[(2.0,"ntransmissions_hyperedge","linear")],
    figs[(0.0,"ntransmissions_hyperedge","sqrt")],
    figs[(1.0,"ntransmissions_hyperedge","sqrt")],
    figs[(2.0,"ntransmissions_hyperedge","sqrt")],
    layout = l,
    size=(2000,800)
)

#remove cbar 
for i in [1,2,4,5]
    plot!(transmission_fig[i],cbar=false)
end
# remove xlabel from first row 
for i in [1,2,3]
    plot!(transmission_fig[i],xlabel="")
end
# remove ylabels 
for i in [2,3,5,6]
    plot!(transmission_fig[i],ylabel="")
end
# adjust margins 
plot!(transmission_fig[4],bottom_margin=8mm,left_margin=14mm)
plot!(transmission_fig[1],top_margin=5mm)
plot!(transmission_fig, plot_title="Transmissions By Hyperedge Size", plot_titlefontsize=20,
        dpi=2000)

savefig(transmission_fig,joinpath(figure_path,"hyperedge_transmissions.png"))
savefig(transmission_fig,joinpath(figure_path,"hyperedge_transmissions.pdf"))



### infections 
l = @layout [grid(2, 3, widths=[0.31, 0.31, 0.38])]
transmission_fig = plot(
    figs[(0.0,"ninfected_hyperedge","linear")],
    figs[(1.0,"ninfected_hyperedge","linear")],
    figs[(2.0,"ninfected_hyperedge","linear")],
    figs[(0.0,"ninfected_hyperedge","sqrt")],
    figs[(1.0,"ninfected_hyperedge","sqrt")],
    figs[(2.0,"ninfected_hyperedge","sqrt")],
    layout = l,
    size=(2000,800)
)

#remove cbar 
for i in [1,2,4,5]
    plot!(transmission_fig[i],cbar=false)
end
# remove xlabel from first row 
for i in [1,2,3]
    plot!(transmission_fig[i],xlabel="")
end
# remove ylabels 
for i in [2,3,5,6]
    plot!(transmission_fig[i],ylabel="")
end
# adjust margins 
plot!(transmission_fig[4],bottom_margin=8mm,left_margin=14mm)
plot!(transmission_fig[1],top_margin=5mm)
plot!(transmission_fig, plot_title="Hyperedge Infected Size", plot_titlefontsize=20,
        dpi=2000)

savefig(transmission_fig,joinpath(figure_path,"hyperedge_infected.png"))
savefig(transmission_fig,joinpath(figure_path,"hyperedge_infected.pdf"))

## Individual Trajectories 
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

f = rolling_hyperedge_transmissions(data)
title!(f,fname)
plot!(f,dpi=1000)
savefig(f,joinpath(figure_path,"rolling_hyperedge_transmissions-$(splitext(fname)[1]).png"))
savefig(f,joinpath(figure_path,"rolling_hyperedge_transmissions-$(splitext(fname)[1]).pdf"))