#### TOTAL INFECTIONS DEMO ####
# for standalone running. loads in aggregated_data dictionaries 
# include("epidemic-figures-utils.jl")

gname = "spatial-hypergraph-50000-2-0.0-1-newalpha.jsonl"
data = aggregated_data[gname]

function stochastic_evolution(data,beta = 9e-1)
    # group data by starting condition, beta, and hyper_beta_func
    filtered_data = hcat(data["beta"], data["hyper_beta_func"], data["initial_num_infected_nodes"])
    grouped_inds = group_by_function(x -> (x[1], x[2], x[3]), eachrow(filtered_data))
    if length(unique(data["delta"]))>1 || length(unique(data["gamma"]))>1
        @warn("mulitple epidemic parameter sets detected. plots likely to exhibit different steady states with respect to beta")
    end

    # put the plot together 
    static_keys = sort(collect(keys(grouped_inds)))
    linear_keys = [k for k in static_keys if k[2] == "linear"]
    sqrt_keys = [k for k in static_keys if k[2] == "sqrt"]

    beta_vals = unique(first.(linear_keys))
    beta_map = Dict(map(x -> (x[2], x[1]), enumerate(sort(beta_vals))))

    # zoom in on beta and a normalization 
    inds = (first.(static_keys).==beta) .& ([x[2]=="linear" for x in static_keys])

    f = Plots.plot(size=(550,400))
    for (i,key) in enumerate(static_keys[inds])
        group = grouped_inds[key]
        Plots.plot!(f,Float64.(data["infections"][group,1:2000]'),c=i,leg=false)
    end
    # plot a square around the trailing infections
    # get all inds for beta and linear normalization 
    group = reduce(union,[grouped_inds[key] for key in static_keys[inds]])
    mean_vals = mean(data["infections"][group,1000:2000],dims=2)
    mean_val = mean(mean_vals)

    offset = 2500
    lw = 5
    Plots.plot!(f,[1000;1000],[mean_val-offset;mean_val+offset], linewidth=lw, c=:black)
    Plots.plot!(f,[2000;2000],[mean_val-offset;mean_val+offset], linewidth=lw, c=:black)
    Plots.plot!(f,[1000;2000],[mean_val+offset;mean_val+offset], linewidth=lw, c=:black)
    Plots.plot!(f,[1000;2000],[mean_val-offset;mean_val-offset], linewidth=lw, c=:black)
    Plots.plot!(f,
        xlabel="Time Step",
        ylabel="Total Infections",
        title = L"\beta = %$beta",
        framestyle=:box,
        thickness_scaling=1.2)
    
    # plot histogram of mean values 
    @show size(mean_vals)
    a,b =extrema(mean_vals)
    println("range for bin: $(b-a)")
    newf = Plots.histogram(vec(mean_vals),nbins=25,
                    title=L"\beta = %$beta",
                    ylabel="Frequency",
                    xlabel="Average Trailing\nTotal Infections",
                    leg=false,
                    size=(500,400))

    return f,newf  
end
#TODO change function call to data after centralizing data loading 
f,newf = stochastic_evolution(data,1e-1)
Plots.plot!(f,
    xscale=:log10,
    # yscale=:log10
    )
display(f)
Plots.plot!(newf,
                xlims=(22300,22500),
                right_margin=5Plots.mm,
                ylims = (0,12),
                framestyle=:box,
                thickness_scaling=1.2,
                margins = 0Plots.mm,
                bottom_margin=-3Plots.mm,
)
display(newf)
Plots.savefig(f,"data/output/figures/final/total-infs-explanation-1.pdf")
Plots.savefig(newf,"data/output/figures/final/total-infs-explanation-2.pdf")
