include("epidemic-figures-utils.jl")

function total_infections_heatmap(agg_data,gname_key="-50000-2-")
    graph_names = collect(keys(aggregated_data))
    graph_names = filter(x->occursin(gname_key,x),graph_names)
    # sort by value of alpha 
    a_vals = map(x->parse(Float64,split(x,"-")[5]),graph_names)
    p = sortperm(a_vals)
    graph_names = graph_names[p]

    b_vals = unique(agg_data[graph_names[1]]["beta"])
    
    # build data for plot
    linear_data = zeros(Float64,length(b_vals),length(ALPHA_VALS))
    pairwise_data = zeros(Float64,length(b_vals),length(ALPHA_VALS))
    sqrt_data = zeros(Float64,length(b_vals),length(ALPHA_VALS))
    
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

beta = 2e-3
f = Plots.scatter(ALPHA_VALS, pairwise_data[BETA_TO_IND[beta],:],
                    markerstrokewidth = 0,
                    markersize=6,
                    yscale=:log10, 
                    leg=false,
                    xlabel = L"\alpha",
                    ylabel="Total Infections",
                    title = L"\beta=%$beta, g(m) = 1", 
                    ylims = (1,25000),
                    thickness_scaling = 1.3,
                )
Plots.savefig(f,"data/output/figures/final/total-infections-50000-2-beta_2e-3-pairwise.pdf")


# second figure 
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

function get_beta_transmissions_data(agg_data,graph_name="-50000-2-",beta=2e-2,gamma=5e-2,delta=1e-2,hyper_beta_func="linear")
    # get related keys 
    graph_names = collect(keys(agg_data))
    graph_names = filter(x->occursin(graph_name,x),graph_names)
    # sort by value of alpha 
    a_vals = map(x->parse(Float64,split(x,"-")[5]),graph_names)
    p = sortperm(a_vals)
    graph_names = graph_names[p]

    # get maximum h_size over all keys and initialize plotting data 
    max_hsize = maximum(map(x->size(agg_data[x]["ntransmissions_hyperedge"],2),graph_names))
    pdata = zeros(Float64,max_hsize,length(ALPHA_VALS))
    # populate data 
    for (col,gname) in enumerate(graph_names)
        data = agg_data[gname]
        # group epidemics and pick out our group
        filtered_data = hcat(data["beta"], data["hyper_beta_func"],data["gamma"],data["delta"])
        grouped_inds = group_by_function(x -> tuple(x...), eachrow(filtered_data))
        inds = grouped_inds[(beta,hyper_beta_func,gamma,delta)]
        # look at infections 
        avg_infs = mean(data["infections"][inds,end-1000:end])
        # populate data 
        target = vec(mean(data["ntransmissions_hyperedge"][inds,:],dims=1))
        for (row,val) in enumerate(target)
            pdata[row,col] = val
        end
    end
    return pdata
end

hyper_beta = "pairwise"
beta = 2e-3
gamma = 5e-2
delta = 5e-2

pdata = get_beta_transmissions_data(aggregated_data,
                                "-50000-2-",
                                beta,
                                gamma,
                                delta,
                                hyper_beta,
                            )
pdata = mapslices(x->x./sum(x), pdata,dims=1)
f = _custom_heatmap(pdata,ALPHA_VALS,false)
Plots.plot!(f,
                xlabel=L"\alpha",
                ylabel="Hyperedge Size",
                title="Normalized Transmissions\nβ=$beta, g(m)=1",
                # clims=(-2.8,0),
                right_margin=3Plots.mm,
                top_margin=3Plots.mm,
                bottom_margin=-1Plots.mm,
                thickness_scaling=1.2)

# Plots.savefig(f,"data/output/figures/final/transmissions-2e-3.pdf")
