using StatsPlots
using Plots
using CairoMakie
using Measures
using StatsBase
using JSON 
using JSON3
using ProgressMeter
import LinearAlgebra: norm
using LaTeXStrings
using DelimitedFiles

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
aggregated_data = load_aggregated_data("aggregated-sirs-output.json")
aggregated_data1 = load_aggregated_data("aggregated-sirs-output-scratch-v1.json")
aggregated_data2 = load_aggregated_data("aggregated-sirs-output-scratch-v2.json")
aggregated_data3 = load_aggregated_data("aggregated-sirs-output-scratch-v3.json")

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

d0
g0

function get_beta_transmissions_data(agg_data,graph_name="-50000-2-",beta=2e-2,gamma=5e-2,delta=1e-2,hyper_beta_func="linear";avg_inf_th=500,coverage=false)
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
        filtered_data = hcat(data["beta"], data["hyper_beta_func"],data["gamma"],data["delta"])
        grouped_inds = group_by_function(x -> tuple(x...), eachrow(filtered_data))
        inds = grouped_inds[(beta,hyper_beta_func,gamma,delta)]
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
        # if coverage
        #     node_coverage = hinfo_dict[gname]["node_coverage"]
        #     node_coverage = vcat(node_coverage,zeros(Float64,size(pdata,1)-length(node_coverage)))
        #     pdata[:,col] ./= node_coverage 
        # end
    end
    return pdata
end

b = 1e-1
g = 5e-2 
hyper_beta = "linear"

d = 5e-2
tmp = new_infections_data(aggregated_data,"-50000-2",b,g,d,hyper_beta)
qinfo = mapslices(x->quantile(x,[0,0.5,1]),tmp,dims=1)'
fig = plot_with_ribbons(qinfo)
Plots.plot!(fig,title=L"\delta=%$d",
                ylabel="Average Total Infections",
                xlabel=L"\alpha")

d = 1e-2
tmp = new_infections_data(aggregated_data2,"-50000-2-",b,g,d,hyper_beta)
qinfo = mapslices(x->quantile(x,[0,0.5,1]),tmp,dims=1)'
fig1 = plot_with_ribbons(qinfo)
Plots.plot!(fig1,title=L"\delta=%$d",
                ylabel="Average Total Infections",
                xlabel=L"\alpha") 
# Plots.savefig(fig1,"data/output/figures/final/hyperedge-effects-p2-1.pdf")               
pdata1 = get_beta_transmissions_data(aggregated_data2,
                                "-50000-2-",
                                b,
                                g,
                                d,
                                hyper_beta,
                                avg_inf_th=0,
                                coverage = false,
                            )
writedlm("data-delta_1e-2.csv", pdata1, ',')
# To load:
# loaded_matrix = readdlm("data-delta_1e-2.csv", ',')
# pdata1[2,:] .= 0
f1 = _custom_heatmap(pdata1,alpha_vals)
Plots.plot!(f1,
                xlabel=L"\alpha",
                ylabel="Hyperedge Size",
                title=L"\delta=%$d",
                clims=(-5.0,0),
                right_margin=3Plots.mm,
                top_margin=3Plots.mm)


d = 2e-2
tmp = new_infections_data(aggregated_data3,"-50000-2-",b,g,d,hyper_beta)
qinfo = mapslices(x->quantile(x,[0,0.5,1]),tmp,dims=1)'
fig2 = plot_with_ribbons(qinfo)
Plots.plot!(fig2,title=L"\delta=%$d",
                ylabel="Average Total Infections",
                xlabel=L"\alpha")
# Plots.savefig(fig2,"data/output/figures/final/hyperedge-effects-p2-2.pdf")               
pdata2 = get_beta_transmissions_data(aggregated_data3,
                                "-50000-2-",
                                b,
                                g,
                                d,
                                hyper_beta,
                                avg_inf_th=0,
                                coverage = false,
                            )
writedlm("data-delta_2e-2.csv", pdata2, ',')
# pdata2[2,:] .= 0
f2 = _custom_heatmap(pdata2,alpha_vals)
Plots.plot!(f2,
                xlabel=L"\alpha",
                ylabel="Hyperedge Size",
                title=L"\delta=%$d",
                clims=(-5.0,0),
                right_margin=3Plots.mm,
                top_margin=3Plots.mm)


d = 3e-2
tmp = new_infections_data(aggregated_data3,"-50000-2-",b,g,d,hyper_beta)
qinfo = mapslices(x->quantile(x,[0,0.5,1]),tmp,dims=1)'
fig3 = plot_with_ribbons(qinfo)
Plots.plot!(fig3,title=L"\delta=%$d",
                ylabel="Average Total Infections",
                xlabel=L"\alpha")
# Plots.savefig(fig3,"data/output/figures/final/hyperedge-effects-p2-3.pdf")               
pdata3 = get_beta_transmissions_data(aggregated_data3,
                                "-50000-2-",
                                b,
                                g,
                                d,
                                hyper_beta,
                                avg_inf_th=0,
                                coverage = false,
                            )
writedlm("data-delta_3e-2.csv", pdata3, ',')
# pdata3[2,:] .= 0
f3 = _custom_heatmap(pdata3,alpha_vals)
Plots.plot!(f3,
                xlabel=L"\alpha",
                ylabel="Hyperedge Size",
                title=L"\delta=%$d",
                clims=(-5.0,0),
                right_margin=3Plots.mm,
                top_margin=3Plots.mm)




extrema(pdata1)
extrema(pdata2)
extrema(pdata3)
Plots.heatmap(abs.(pdata2.-pdata1),
    yscale=:log10,
    clims=(0,150))


Plots.heatmap(sqrt.(abs.(pdata3.-pdata1)),
    yscale=:log10,
    # clims=(0,100),
    c=colormap("RdBu"))

_custom_heatmap(abs.(pdata3.-pdata2),alpha_vals,true)


pdata1
pdata2.-pdata1
pdata3.-pdata1
pr3 = mapslices(x->x./sum(x),pdata3,dims=1)
pr2 = mapslices(x->x./sum(x),pdata2,dims=1)
pr1 = mapslices(x->x./sum(x),pdata1,dims=1)


f = _custom_heatmap(abs.(pr2.-pr1),alpha_vals,false)
Plots.plot!(f,clims=(-6,-1.3), title=L"(\delta_{2e-2})-(\delta_{1e-2})")
f = _custom_heatmap(abs.(pr3.-pr1),alpha_vals,false)
Plots.plot!(f,clims=(-6,-1.3), title=L"(\delta_{3e-2})-(\delta_{1e-2})")


f = _custom_heatmap(abs.(pdata2.-pdata1),alpha_vals,false)
Plots.plot!(f,clims=(-4,2.5), title=L"(\delta_{2e-2})-(\delta_{1e-2})")
f = _custom_heatmap(abs.(pdata3.-pdata1),alpha_vals,false)
Plots.plot!(f,clims=(-4,2.5), title=L"(\delta_{3e-2})-(\delta_{1e-2})")

pdata



mylayout = [Plots.grid(3, 2, widths=[0.45, 0.55])]
Plots.plot(fig1,f1,fig2,f2,fig3,f3,
                layout = mylayout,
                size=(1600,1600),
                thickness_scaling=1.2,
            # link=:y,
            # yscale=:log10,
            # ylims=(6000,17000)
)
f1
f2
f3
# plot 1: can we see hysteresis as we vary beta for a single graph 
# plot binned rolling average by starting condition vs beta as scatter plot
# if there is no bistability (not in that parameter regime), then
# we'd expect to see that points corresponding to the same beta should be very close to one another
# (the initial condition doesn't matter)
# otherwise, we'd expect to see meaningfully different values for the same beta 
using Plots 
function _custom_heatmap(pdata,xvals = -1,max_normalize=true)
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
         size = (400, 300))
end

graph_names = collect(keys(aggregated_data))
graph_names = filter(x->occursin("-50000-2-",x),graph_names)
graph_names = sort(graph_names,by=x->parse(Float64,split(x,"-")[5]))

data = aggregated_data[graph_names[1]]
unique(data["beta"])
unique(data["gamma"])
unique(data["delta"])

g,d = 5e-2,1e-2
hyper_beta = "linear"
for b in unique(data["beta"])
    tmp = new_infections_data(aggregated_data3,"-50000-2",b,g,d,hyper_beta)
    qinfo = mapslices(x->quantile(x,[0,0.5,1]),tmp,dims=1)'
    fig = plot_with_ribbons(qinfo)
    Plots.plot!(fig,title="beta=$b,delta=$d, g(m)=$hyper_beta",
                    xlabel="alpha")
    display(fig)
end


g,d = 5e-2, 1e-2
figs = []
for b in [1e-2,2e-2,3e-2, 4e-2, 5e-2,1e-1,9e-1]
    tmp = new_infections_data(aggregated_data3,"-50000-2-",b,g,d,"linear")
    qinfo = mapslices(x->quantile(x,[0,0.5,1]),tmp,dims=1)'

    lower = qinfo[:, 1]
    middle = qinfo[:, 2]
    upper = qinfo[:, 3]

    f = Plots.plot(alpha_vals, middle, 
        ribbon = (middle - lower, upper - middle),
        fillalpha = 0.3,
        linewidth = 2,
        color = :blue,
        xlabel = L"\alpha",
        ylabel = "Average Total Infections",
        title = L"\beta=%$b",
        legend = false,
        size = (400, 300))
    # display(f)
    push!(figs,f)
end 
f = Plots.plot(figs[[1,3,6,7]]...,
        layout=(1,4),
        size=(1400,300),
        # link=:y,
)

for (ind,f) in enumerate(figs[[1,3,6,7]])
    Plots.savefig(f,"data/output/figures/final/hyperedge-effects-p1-$ind.pdf")
end


# Plots.plot!(yscale=:log10,ylims=(10,25000))
    
pdata1
graph_names = sort(graph_names,by=x->parse(Float64,split(x,"-")[5]))

graph_names[1]

# specify epidemic params and show a single trace (trajectory with bands)
data = aggregated_data[graph_names[1]]
filtered_data = hcat(data["beta"],data["gamma"],data["delta"],data["hyper_beta_func"],data["initial_num_infected_nodes"])
grouped_inds = group_by_function(x->tuple(x...),eachrow(filtered_data))

first(keys(grouped_inds))
group = (1e-1,5e-2,4e-2,"linear",5000)
inds = grouped_inds[group]

Array(Float64.(data["infections"][inds,:]'))
f = Plots.plot(size=(300,300))
Plots.plot!(Array(Float64.(data["infections"][inds,:]')),
        xscale=:log10,
        yscale=:log10,
        leg=false,
        c=:blue
)
f = Plots.plot(size=(300,300))
Plots.plot!(Array(Float64.(data["infections"][inds,end-1000:end]')),
        # xscale=:log10,
        # yscale=:log10,
        leg=false,
        c=:blue
)
mean(Array(Float64.(data["infections"][inds,end-1000:end]')),dims=1)



# smaller graph 
graph_names = collect(keys(aggregated_data))
data = aggregated_data[graph_names[1]]
unique(data["beta"])
unique(data["gamma"])
unique(data["delta"])


g,d = 5e-2, 4e-2
figs = []
for b in [1e-2, 2e-2, 3e-2, 4e-2, 5e-2, 1e-1, 9e-1]
    tmp = new_infections_data(aggregated_data,"-50000-2-",b,g,d,"linear")
    qinfo = mapslices(x->quantile(x,[0,0.5,1]),tmp,dims=1)'

    lower = qinfo[:, 1]
    middle = qinfo[:, 2]
    upper = qinfo[:, 3]

    f = Plots.plot(alpha_vals, middle, 
        ribbon = (middle - lower, upper - middle),
        fillalpha = 0.3,
        linewidth = 2,
        color = :blue,
        xlabel = "alpha",
        ylabel = "Average Infections",
        title = "beta=$b, delta =$d\ngamma=$g",
        legend = false,
        size = (400, 300))
    # display(f)
    push!(figs,f)
end 

f = Plots.plot(figs[[3,4,5,6]]...,
                layout=(1,4),
                size=(1400,300),
)

display(f)






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
        i0 = unique(data["initial_num_infected_nodes"])
        # group data 
        filtered_data = hcat(data["beta"], data["hyper_beta_func"], data["initial_num_infected_nodes"])
        grouped_inds = group_by_function(x -> (x[1], x[2], x[3]), eachrow(filtered_data))
        grouped_averages = Dict(k => mean(data["infections"][grouped_inds[k],:]) for k in keys(grouped_inds))
        for (row,beta) in enumerate(beta_vals)
            linear_data[row,col] = _compute_range([grouped_averages[(beta,"linear",i)] for i in i0])
            sqrt_data[row,col] = _compute_range([grouped_averages[(beta,"sqrt",i)] for i in i0])
            pairwise_data[row,col] = _compute_range([grouped_averages[(beta,"pairwise",i)] for i in i0])
        end
    end    
    return linear_data, sqrt_data, pairwise_data    
end



graph_names = collect(keys(aggregated_data3))
graph_names = filter(x->occursin("-50000-2-",x),graph_names)
graph_names = sort(graph_names,by=x->parse(Float64,split(x,"-")[5]))

#### TOTAL INFECTIONS DEMO ####
tmp_gname = "spatial-hypergraph-50000-2-0.0-1-newalpha.jsonl"
data = load_epidemic_data("data/hysteresis/sirs/",tmp_gname,excluded=["ninfected_hyperedge","ntransmissions_hyperedge"],truncate=false)


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



pdata1




_, hist1, ybins1 = _custom_heatmap(pdata1,alpha_vals,false,true)
_, hist2, ybins2 = _custom_heatmap(pdata2,alpha_vals,false,true)
_, hist3, ybins3 = _custom_heatmap(pdata3,alpha_vals,false,true)

function max_normalize_col(data)
    return mapslices(x->x./maximum(x),data,dims=1)
end

function sum_normalize_col(data)
    return mapslices(x->x./sum(x),data,dims=1)
end


Plots.heatmap(sum_normalize_col(hist1) ./ sum_normalize_col(hist2))