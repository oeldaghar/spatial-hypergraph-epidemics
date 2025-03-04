using Plots, LaTeXStrings, Colors
using Measures
using JSON 

# load data in 
n = 10000
println("Reading data...")
loaded_data = open("data/output/alpha_func_data-n_$n.json", "r") do file
    JSON.parse(file)
end

newdata = []
push!(newdata,loaded_data["($n, 2, \"alpha_func\")"])
push!(newdata,loaded_data["($n, 2, \"linear\")"])
push!(newdata,loaded_data["($n, 2, \"no_dim\")"])
push!(newdata,loaded_data["($n, 5, \"alpha_func\")"])
push!(newdata,loaded_data["($n, 5, \"linear\")"])
push!(newdata,loaded_data["($n, 5, \"no_dim\")"])
push!(newdata,loaded_data["($n, 10, \"alpha_func\")"])
push!(newdata,loaded_data["($n, 10, \"linear\")"])
push!(newdata,loaded_data["($n, 10, \"no_dim\")"])

newdata = map(x->reduce(hcat,x),newdata)


function quantile_figure(data,alphas)
    quantiles = [0.1,0.25,0.5,0.75,0.9]
    linewidths = [0.5, 1, 2, 1, 0.5]
    colors = [:lightgrey, :grey, :black, :grey, :lightgrey]

    # total edges 
    f = Plots.plot(xlabel=L"\alpha", ylabel="Total hyperedges", 
                framestyle=:box,
                thickness_scaling=1.2,
                guidefontsize=14,
                tickfontsize=12,
                tickdirection=:out)
    for (q, lw, c) in zip(quantiles, linewidths, colors)
        nedges_q = quantile.(eachrow(data), q)
        Plots.plot!(f, alphas, nedges_q, label="", linewidth=lw, color=c,    
            yscale=:log10)
    end
    return f
end

function _remove_tick_labels(f)
    curr_yticks = Plots.yticks(f)
    new_yticks = (curr_yticks[1],["" for i=1:lastindex(curr_yticks[1])])
    Plots.plot!(f,yticks=new_yticks,ylabel="")
end    

figs = map(x->quantile_figure(x,range(0,2,25)),newdata)

# put figure together 
plt = Plots.plot(figs...,layout = (3,3),
                    size=(1500,1400),
                    top_margin=8Measures.mm,
                    link=:all)
# touch up margins and label 
Plots.plot!(plt[2],title = "n=$n, d=2",
        titlefontsize=22)
Plots.plot!(plt[5],title = "n=$n, d=5",
        titlefontsize=22)
Plots.plot!(plt[8],title = "n=$n, d=10",
        titlefontsize=22)
for i=1:9
    if i%3!=1
        _remove_tick_labels(plt[i])
    end
end
Plots.plot!(plt,bottom_margin=2mm)
Plots.plot!(plt[1],left_margin=8mm)
Plots.plot!(plt,top_margin=5mm,dpi = 1000)

Plots.savefig(plt,"data/output/figures/final/alpha-figure.pdf")