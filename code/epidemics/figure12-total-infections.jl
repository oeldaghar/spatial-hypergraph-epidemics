# for standalone running. loads in aggregated_data dictionaries  
# include("epidemic-figures-utils.jl")


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
        grouped_averages = Dict(k => mean(data["infections"][grouped_inds[k],end-1000:end]) for k in keys(grouped_inds))
        for (row,beta) in enumerate(b_vals)
            linear_data[row,col] = grouped_averages[(beta,"linear")]
            sqrt_data[row,col] = grouped_averages[(beta,"sqrt")]
            pairwise_data[row,col] = grouped_averages[(beta,"pairwise")]
        end
    end    
    return linear_data, sqrt_data, pairwise_data
end

linear_data, sqrt_data, pairwise_data = total_infections_heatmap(aggregated_data,"-50000-2-")

# for colorbar limits 
@show extrema(log10.(linear_data))
@show extrema(log10.(sqrt_data))
@show extrema(log10.(pairwise_data))

# first figure.. heatmaps 
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


# second figure.. trajectories
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
