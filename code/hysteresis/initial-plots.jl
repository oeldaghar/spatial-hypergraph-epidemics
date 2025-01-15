# initial figures for hysteresis effects 

using Plots
using Measures
using StatsBase
using JSON 

data_dir = "data/hysteresis/sirs/"

gname = "spatial-hypergraph-5000-2"
fnames = filter(x->occursin(gname,x),readdir(data_dir))

alphas = map(x->parse(Float64,split(x,"-")[5]),fnames)
perm = sortperm(alphas)
alphas = alphas[perm]
fnames = fnames[perm]

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

# plot 1: can we see hysteresis as we vary beta for a single graph 
fnames[1]
fnames[8]
fnames[end]
fname = fnames[1]
# plot binned rolling average by starting condition vs beta as scatter plot
# if there is no bistability (not in that parameter regime), then
# we'd expect to see that points corresponding to the same beta should be very close to one another
# (the initial condition doesn't matter)
# otherwise, we'd expect to see meaningfully different values for the same beta 

# Example usage
function plot_hysteresis(data_dir, fname)
    # alpha_val = Float64(split(fname,"-")[5])
    data = JSON.parsefile(joinpath(data_dir, fname))
    data["beta"] = Float64.(data["beta"])
    data["hyper_beta_func"] = String.(data["hyper_beta_func"])
    data["infections"] = hcat(Vector{Int}.(data["infections"])...)

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
                xlabel = "Pairwise Infection Probability", 
                ylabel = "Binned Average\nRolling Infections", 
                title = "Beta vs Rolling Infections - $title_suffix\n$fname",
                markerstrokewidth = 0,
                xlims=(8e-4,1.05e0))
        plot!(f,xscale=:log10)
    end

    return plot_keys(linear_keys, "Linear"),plot_keys(sqrt_keys, "Sqrt")
end

# Example usage
# data_dir = "data/hysteresis/sirs/"
# fname = "spatial-hypergraph-5000-2-0.0-1.json"
# f1,f2 = plot_hysteresis(data_dir, fname)
# f1

fs = []
for fname in [
            "spatial-hypergraph-5000-2-0.0-1.json",
            "spatial-hypergraph-5000-2-1.0-1.json",
            "spatial-hypergraph-5000-2-2.0-1.json"
        ]
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
f_main

# save figure 
data_path = "data/hysteresis/figures"
if !ispath(data_path)
    mkpath(data_path)
end
savefig(f_main, joinpath(data_path,"hysteresis_plot_v1.pdf"))
plot!(f_main,dpi=1000)
savefig(f_main, joinpath(data_path,"hysteresis_plot_v1.png"))

#zooming in on one piece 
function stochastic_evolution(data_dir,fname)
    data = JSON.parsefile(joinpath(data_dir, fname))
    data["beta"] = Float64.(data["beta"])
    data["hyper_beta_func"] = String.(data["hyper_beta_func"])
    data["infections"] = hcat(Vector{Int}.(data["infections"])...)

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
        title = "alpha=$alpha_val, beta = $beta")
    return f 
end

filtered_fnames = [
            "spatial-hypergraph-5000-2-0.0-1.json",
            "spatial-hypergraph-5000-2-1.0-1.json",
            "spatial-hypergraph-5000-2-2.0-1.json"
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
savefig(main_f, joinpath(data_path,"hysteresis_stochastic_plot_v1.pdf"))
plot!(main_f,dpi=1000)
savefig(main_f, joinpath(data_path,"hysteresis_stochastic_plot_v1.png"))
        