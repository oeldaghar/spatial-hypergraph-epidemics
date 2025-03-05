
# for standalone running. loads in aggregated_data dictionaries  
# include("epidemic-figures-utils.jl")

function hysteresis_data(agg_data,gname_key="-50000-2-")
    graph_names = collect(keys(agg_data))
    graph_names = filter(x->occursin(gname_key,x),graph_names)
    # sort by value of alpha 
    alpha_vals = map(x->parse(Float64,split(x,"-")[5]),graph_names)
    p = sortperm(alpha_vals)
    graph_names = graph_names[p]
    

    # build data for plot
    linear_data = zeros(Float64,length(BETA_VALS),length(ALPHA_VALS))
    pairwise_data = zeros(Float64,length(BETA_VALS),length(ALPHA_VALS))
    sqrt_data = zeros(Float64,length(BETA_VALS),length(ALPHA_VALS))
    
    function _compute_range(itr)
        min_val,max_val = extrema(itr)
        return max_val - min_val
    end

    for (col,gname) in enumerate(graph_names) # loop over alpha
        data = agg_data[gname]
        # group data 
        filtered_data = hcat(data["beta"], data["hyper_beta_func"])
        grouped_inds = group_by_function(x -> tuple(x...), eachrow(filtered_data))
        for (row,beta) in enumerate(BETA_VALS)
            linear_data[row,col] = _compute_range(mean(data["infections"][grouped_inds[(beta,"linear")],end-1000:end],dims=2))
            sqrt_data[row,col] = _compute_range(mean(data["infections"][grouped_inds[(beta,"sqrt")],end-1000:end],dims=2))
            pairwise_data[row,col] = _compute_range(mean(data["infections"][grouped_inds[(beta,"pairwise")],end-1000:end],dims=2))
        end
    end    
    return linear_data, sqrt_data, pairwise_data    
end

linear_hys, sqrt_hys, pairwise_hys = hysteresis_data(aggregated_data,"-50000-2-")

# for setting the maximum on the colorbar 
@show maximum(linear_hys)
@show maximum(sqrt_hys)
@show maximum(pairwise_hys)

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
    f = Plots.heatmap(ALPHA_VALS, BETA_VALS, tmp,
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
