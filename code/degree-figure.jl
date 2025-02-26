using JSON
using Plots
using Measures


function plot_data(trials)
    nedges = reduce(hcat, trials) # alpha as row index and trial as column index
    # compute quantiles over columns
    quantiles = [0.1,0.25,0.5,0.75,0.9]
    linewidths = [0.5, 1, 2, 1, 0.5]
    colors = [:lightgrey, :grey, :black, :grey, :lightgrey]
    p = Plots.plot(xlabel="α", ylabel="Weighted Projected\nEdge Volumne",
        # framestyle=:grid,
        tickfontsize=10,
        yguidefontsize = 15,
        xguidefontsize = 15)
    for (q, lw, c) in zip(quantiles, linewidths, colors)
        nedges_q = quantile.(eachrow(nedges), q)
        Plots.plot!(p, alphas, nedges_q, label="", linewidth=lw, color=c)
    end
    p
end


#load and parse data 
projected_trials_pairwise = JSON.parsefile("data/output/spatial_graph_projected_degree_50000_2_pairwise.json")["projected_trials"]
projected_trials_sqrt = JSON.parsefile("data/output/spatial_graph_projected_degree_50000_2_sqrt.json")["projected_trials"]
projected_trials_linear = JSON.parsefile("data/output/spatial_graph_projected_degree_50000_2_linear.json")["projected_trials"]

projected_trials_pairwise = Vector{Vector{Float64}}(projected_trials_pairwise)
projected_trials_sqrt = Vector{Vector{Float64}}(projected_trials_sqrt)
projected_trials_linear = Vector{Vector{Float64}}(projected_trials_linear)

#standalone plot
p1 = plot_data(projected_trials_pairwise)
Plots.plot!(p1,title="g(m)=1",
    titlefont=font("Helvetica Bold", 14))

p2 = plot_data(projected_trials_sqrt)
Plots.plot!(p2,title="g(m)=sqrt(m)",
    titlefont=font("Helvetica Bold", 14))

p3 = plot_data(projected_trials_linear)
Plots.plot!(p3,title="g(m)=m",
    titlefont=font("Helvetica Bold", 14))

#### COMBINING PLOTS
figs = Plots.plot(p3,p2,p1,
            layout=(1,3),size=(1500,400),
            bottom_margin=10Measures.mm,
            # link=:y,
)
#visual spacing
Plots.plot!(figs[1],left_margin=14Measures.mm)
Plots.plot!(figs[2],left_margin=8Measures.mm)
Plots.plot!(figs[3],left_margin=8Measures.mm)
#increase resolution
Plots.plot!(figs,titlefontsize=20,tickfontsize=15,guidefontsize=16,dpi=300)

# Plots.savefig(figs,"data/output/figures/final/projected-degree.pdf")
# Plots.savefig(figs,"data/output/figures/poster/projected-degree.svg")

# #link yaxes for a version 2 of the plot
# Plots.plot!(figs,link=:y)
# ticks = yticks(figs[1])
# Plots.plot!(figs[2], ylabel="",framestyle=:grid,
#             yticks=(ticks[1],fill("",length(ticks[2]))))
# Plots.plot!(figs[3], ylabel="",framestyle=:grid,
#             yticks=(ticks[1],fill("",length(ticks[2]))))
# Plots.plot!(figs,background_color_inside=:gray98,top_margin=3Measures.mm)
# Plots.savefig(figs,"data/output/figures/final/projected-degree-v2.pdf")
# Plots.savefig(figs,"data/output/figures/poster/projected-degree-v2.svg")