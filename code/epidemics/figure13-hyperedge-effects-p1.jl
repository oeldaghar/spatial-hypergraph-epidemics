# for standalone running. loads in aggregated_data dictionaries  
# include("epidemic-figures-utils.jl")
# aggregated_data3 = load_aggregated_data("aggregated_data","aggregated-sirs-output-scratch-v3.json")

function new_infections_data(agg_data::Dict,
                                gname_key::String,
                                beta::Float64,
                                gamma::Float64,
                                delta::Float64,
                                hyper_beta_func::String)
    graph_names = collect(keys(agg_data))
    graph_names = filter(x->occursin(gname_key,x),graph_names)
    # sort by value of alpha 
    a_vals = map(x->parse(Float64,split(x,"-")[5]),graph_names)
    p = sortperm(a_vals)
    graph_names = graph_names[p]

    inf_data = Vector()
    for gname in graph_names
        data = agg_data[gname]
        filtered_data = hcat(data["beta"], data["gamma"], data["delta"], data["hyper_beta_func"])
        grouped_inds = group_by_function(x->tuple(x...),eachrow(filtered_data))

        inds = grouped_inds[(beta,gamma,delta,hyper_beta_func)]
        push!(inf_data,rolling_average.(eachrow(data["infections"][inds,:]),1000))
    end
    return reduce(hcat,inf_data)
end


gamma,delta = 5e-2, 1e-2
hyper_beta_func = "linear"
figs = []
for (ind,b) in enumerate([1e-2,3e-2,1e-1,9e-1])
    tmp = new_infections_data(aggregated_data3,"-50000-2-",b,gamma,delta,hyper_beta_func)
    qinfo = mapslices(x->quantile(x,[0,0.5,1]),tmp,dims=1)'

    lower = qinfo[:, 1]
    middle = qinfo[:, 2]
    upper = qinfo[:, 3]

    f = Plots.plot(ALPHA_VALS, middle, 
        ribbon = (middle - lower, upper - middle),
        fillalpha = 0.3,
        linewidth = 2,
        color = :blue,
        xlabel = L"\alpha",
        ylabel = "Average Total Infections",
        title = L"\beta=%$b",
        legend = false,
        size = (400, 300))
    display(f)
    Plots.savefig(f,"data/output/figures/final/hyperedge-effects-p1-$ind.pdf")
    push!(figs,f)
end 

# side by side 
# f = Plots.plot(figs...,
#         layout=(1,4),
#         size=(1400,300),
#         # link=:y,
# )