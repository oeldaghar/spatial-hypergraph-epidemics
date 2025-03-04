# Plotting 
using Plots, LaTeXStrings, Colors
using LazySets
using Measures
using JSON 


data = JSON.parsefile("data/output/hyper_ppr_data.json")
X = data["X"]
X = map(x->[x[1]; x[2]],X)
X = reduce(hcat,X)

graphs = data["graphs"]
graphs = Vector{Vector{Int}}.(graphs)

ppr_soln = data["ppr_soln"]
ppr_soln = Vector{Float64}.(ppr_soln)

alphas = Float64.(data["alphas"])

function plot_hypergraph(X,hedges,nodes)
    plt = Plots.plot(leg=false, size=(300,300))
    for edge in hedges
        if length(edge) > 2 
            local b = [Ball2(X[:, v], 0.0005) for v in edge]
            local c = ConvexHullArray(b)
            Plots.plot!(plt, c, 1e-3, alpha=0.05, c=:blue)
        else
            Plots.plot!(plt, X[1,edge], X[2,edge], c=:blue,alpha=0.35)
        end
    end
    return plt
end

function make_ppr_plot(hedges,ppr)
    new_inds = findall(ppr.>0)

    newx = log10.(1e-16 .+ ppr)
    markersizes = 3 .*ones(Float64,lastindex(newx))
    markersizes[ppr.>0] .= 8

    plt = plot_hypergraph(X,hedges,1:size(X,2))
    Plots.scatter!(plt,X[1,:],X[2,:],marker_z=newx,
            markerstrokewidth=0.5,
            markersize=markersizes,
            xlabel = "Spatial Coordinate 1",
            ylabel = "Spatial Coordinate 2",
            framestyle=:box,
            )
    Plots.scatter!(plt,X[1,new_inds],X[2,new_inds],marker_z=newx[new_inds],
        markerstrokewidth=0.5,
        markersize=markersizes[new_inds],
        xlabel = "Spatial Coordinate 1",
        ylabel = "Spatial Coordinate 2",
        framestyle=:box,
        )
    return plt 
end

figs = map(x->make_ppr_plot(x...),zip(graphs,ppr_soln))

for (ind,f) in enumerate(figs)
    Plots.plot!(f,xlims=(0.62,1.01), ylims=(0.62,1.01),
            title=L"\alpha = %$(alphas[ind])",
            titlefontsize=18,
            bottom_margin=0Plots.mm,
            top_margin = 0Plots.mm,
            left_margin=0Plots.mm)
    Plots.savefig(f,"data/output/figures/final/ppr-alpha_$(alphas[ind]).pdf")    
    display(f)
end

# construct larger plot 
ppr_soln = reduce(hcat,ppr_soln)
f = Plots.plot(size=(500,400))
Plots.plot!(f,[-16;0],[-16;0],
            linewidth = 2,
            color = :black,
            label="",
            legendfontsize=12)
colors = [:grey, :blue, :green]
for col=1:lastindex(ppr_soln,2)
    Plots.scatter!(f,log10.(ppr_soln[:,1]),log10.(ppr_soln[:,col]),
            markerstrokewidth=0,
            markersize=5,
            c=colors[col],
            label=L"\alpha=%$(alphas[col])"
            )
end
minval = minimum(x-> !isinf(x) ? x : 1,log10.(ppr_soln))
Plots.plot!(f,
            xlims=(minval-1e-1,0),
            ylims=(minval-1e-1,0),
            xlabel = "Pairwise PPR Solution",
            ylabel = "Hypergraph PPR Solutions",
            title="PPR Solutions",
)
# manually update ticks 
new_xticks = xticks(f)[1]
new_xticks = (new_xticks[1], map(x->"10^{$x}",new_xticks[2]))
new_yticks = yticks(f)[1]
new_yticks = (new_yticks[1], map(x->"10^{$x}",new_yticks[2]))
Plots.plot!(f, xticks=new_xticks, yticks = new_yticks, framestyle=:box,
        thickness_scaling=1.1)
Plots.savefig(f,"data/output/figures/final/ppr-solutions.pdf")