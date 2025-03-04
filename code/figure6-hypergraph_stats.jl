# Plotting 
using Plots, LaTeXStrings, Colors
using Measures

using JSON 

# load in the data 
n = 10000
json_data = JSON.parsefile("data/output/hypergraph_stats-n_$n.json")

# Function to reconstruct a row
function reconstruct_row(n, k)
    edge_cnts = json_data["($n, $k, \"edge_cnts\")"]
    nedges = json_data["($n, $k, \"nedges\")"]
    ntris = json_data["($n, $k, \"ntris\")"]
    
    # Convert arrays to matrices
    edge_cnts_matrix = reduce(hcat, edge_cnts)
    nedges_matrix = reduce(hcat, nedges)
    ntris_matrix = reduce(hcat, ntris)
    
    return [edge_cnts_matrix, nedges_matrix, ntris_matrix]
end

# Reconstruct rows
row1 = reconstruct_row(n, 2)
row2 = reconstruct_row(n, 5)
row3 = reconstruct_row(n, 10)

# Combine rows into data
data = [row1, row2, row3]


function _custom_heatmap(pdata)
    # function to bin rows of the heatmap weight matrix
    function _bin_indices(bins::Vector, pdata::Vector)
        bin_dict = Dict{Float64, Vector{Int}}()
        for i in 1:length(bins)-1
            bin_dict[bins[i]] = findall(x -> bins[i] <= x < bins[i+1], pdata)
        end
        bin_dict[bins[end]] = findall(x -> x >= bins[end], pdata)
        return bin_dict
    end
    
    function _log_bin_ydata(pdata,ybins)
        binned_inds = _bin_indices(ybins,collect(1:lastindex(pdata,1)))
        new_mat = zeros(Float64,lastindex(ybins)-1,lastindex(pdata,2))
        for (newrow,key) in enumerate(ybins[1:end-1])
            old_rows = binned_inds[key]
            new_mat[newrow,:] = sum(pdata[old_rows,:],dims=1)
        end
        return new_mat
    end
    
    max_hsize = lastindex(pdata,1)
    ybins = (10.0).^(range(1,log10(max_hsize+10),15))
    ybins = vcat(1:9,ybins)
    ybins = sort(unique(ybins))

    new_mat = _log_bin_ydata(pdata,ybins)

    # xrange, yrange, data_matrix, yaxis scale
    f = Plots.heatmap(1:lastindex(new_mat,2),ybins[1:end-1],log10.(new_mat),
                yscale=:log10,
                color=:viridis,
                # clims=(-5,0),
                )
    return f#,new_mat,ybins 
end

function make_fig(data)
    # handle heatmaps
    # make same size 
    max_hyperedge_size = maximum(map(x->size(x[1],1),data))
    for row in data
        heatmap_data = row[1]
        rows_to_pad = max_hyperedge_size-size(heatmap_data,1)
        row[1] = vcat(heatmap_data,zeros(rows_to_pad,size(heatmap_data,2)))
    end

    # make individual figures 
    # heatmaps 
    col1_figs = []
    for pdata in first.(data)
        f = _custom_heatmap(pdata)
        Plots.plot!(f,
                xlabel=L"\alpha",
                ylabel="Hyperedge size",
                colorbar_titlefontsize = 14,
                xticks=([1,7,13,19,25], [0.0, 0.5, 1.0, 1.5, 2.0]),
                yscale=:log10,
                color=:viridis,
                # clims=(-6.0,0),
                framestyle=:box,
                thickness_scaling=1.2,
                guidefontsize=14,
                tickfontsize=12)
        push!(col1_figs,f)
    end

    # col2 and col3 
    quantiles = [0.1,0.25,0.5,0.75,0.9]
    linewidths = [0.5, 1, 2, 1, 0.5]
    colors = [:lightgrey, :grey, :black, :grey, :lightgrey]

    # total edges 
    col2_figs = []
    for pdata in map(x->x[2],data)
        f = Plots.plot(xlabel=L"\alpha", ylabel="Total hyperedges", 
                    framestyle=:box,
                    thickness_scaling=1.2,
                    guidefontsize=14,
                    tickfontsize=12,
                    tickdirection=:out,
                    )
        for (q, lw, c) in zip(quantiles, linewidths, colors)
            nedges_q = quantile.(eachrow(pdata), q)
            Plots.plot!(f, alphas, nedges_q, label="", linewidth=lw, color=c,    
                yscale=:log10)
        end
        push!(col2_figs,f)
    end
    # align ylims 
    col_ylims = Plots.ylims(Plots.plot(col2_figs...,layout=(3,1),link=:all))
    for f in col2_figs
        Plots.plot!(f,ylims=col_ylims)
    end

    col3_figs = []
    for pdata in map(x->x[3],data)
        f = Plots.plot(xlabel=L"\alpha", ylabel="Total Triangles", 
                    framestyle=:box,
                    thickness_scaling=1.2,
                    guidefontsize=14,
                    tickfontsize=12,
                    tickdirection=:out,
                    )
        for (q, lw, c) in zip(quantiles, linewidths, colors)
            ntris_q = quantile.(eachrow(pdata), q)
            Plots.plot!(f, alphas, ntris_q, label="", linewidth=lw, color=c,    
                yscale=:log10)
        end
        push!(col3_figs,f)
    end
    # align ylims 
    col_ylims = Plots.ylims(Plots.plot(col3_figs...,layout=(3,1),link=:all))
    for f in col3_figs
        Plots.plot!(f,ylims=col_ylims)
    end

    # put them all together 
    figs =[]
    for tup in zip(col1_figs,col2_figs,col3_figs)
        push!(figs,tup[1])
        push!(figs,tup[2])
        push!(figs,tup[3])
    end

    l = @layout [Plots.grid(3, 3, widths=[0.42, 0.29, 0.29])]
    plt = Plots.plot(figs...,layout=l, 
                margin=0*Plots.mm, size=(1200,1100))
    Plots.plot!(plt,top_margin = 10mm,bottom_margin=-2mm)
    Plots.plot!(plt[2],title="n=$n  d=2",titlefontsize = 20,
                    top_margin=-3Measures.mm)
    Plots.plot!(plt[5],title="n=$n  d=5",titlefontsize = 20,
                    top_margin=-3Measures.mm)
    Plots.plot!(plt[8],title="n=$n  d=10",titlefontsize = 20,
                    top_margin=-3Measures.mm)
    return plt
end

plt = make_fig(data)
Plots.savefig(plt,"data/output/figures/final/hypergraph-stats.pdf")