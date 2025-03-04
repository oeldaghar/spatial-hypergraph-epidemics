using JSON
using Plots
using Measures
using LaTeXStrings

function plot_data(trials)
    eigvals = reduce(hcat, trials) # now we have alpha as row index and trial as column index
    p = Plots.plot(xlabel="α", ylabel="Weighted Projected λ₁",
                    tickfontsize=10,
                    yguidefontsize = 15,
                    xguidefontsize = 15)  

    # show the full data as a ribbon
    minl = minimum.(eachrow(eigvals))
    maxl = maximum.(eachrow(eigvals))
    mid = (minl .+ maxl)/2
    upper = maxl .- mid
    lower = mid .- minl
    @show minl, maxl
    alphas = range(0,2,lastindex(trials,1))
    Plots.plot!(p, alphas, (minl .+ maxl)/2, 
                linewidth=0, 
                linealpha=0,
                ribbon=(upper, lower), 
                label="",
                c=Gray(0.8))

    # I want to take the quantiles over the columns... 
    quantiles = [0.1,0.5,0.9]
    linewidths = [1.5, 3, 1.5]
    colors = [ :grey, :black, :grey,]

    for (q, lw, c) in zip(quantiles, linewidths, colors)
        eigvals_q = quantile.(eachrow(eigvals), q)
        Plots.plot!(p, alphas, eigvals_q, label="", linewidth=lw, color=c)
    end

    p
end
  
#load and parse data 
projected_trials_pairwise = JSON.parsefile("spatial_graph_projected_lambda1_50000_2_pairwise.json")["projected_trials"]
projected_trials_sqrt = JSON.parsefile("spatial_graph_projected_lambda1_50000_2_sqrt.json")["projected_trials"]
projected_trials_linear = JSON.parsefile("spatial_graph_projected_lambda1_50000_2_linear.json")["projected_trials"]

projected_trials_pairwise = Vector{Vector{Float64}}(projected_trials_pairwise)
projected_trials_sqrt = Vector{Vector{Float64}}(projected_trials_sqrt)
projected_trials_linear = Vector{Vector{Float64}}(projected_trials_linear)

#standalone plots
p1 = plot_data(projected_trials_pairwise)
# ticks = yticks(p1)[1]
# Plots.plot!(p1, ylabel="",framestyle=:grid,
#             yticks=(ticks[1],fill("",length(ticks[2]))))
Plots.plot!(p1,title="g(m)=1",
    titlefont=font("Helvetica Bold", 14),
    ylims = (7.548591979688312, 879.3020219519976),
    yscale=:log10,
    )
Plots.savefig(p1,"data/output/figures/final/projected-lambda1-pairwise.pdf")

p2 = plot_data(projected_trials_sqrt)
Plots.plot!(p2,title="g(m)=sqrt(m)",
    titlefont=font("Helvetica Bold", 14),
    ylims = (7.548591979688312, 879.3020219519976),
    yscale=:log10)
Plots.savefig(p2,"data/output/figures/final/projected-lambda1-sqrt.pdf")


p3 = plot_data(projected_trials_linear)
Plots.plot!(p3, title="g(m)=m",
            titlefont=font("Helvetica Bold", 14),
            ylims = (7.548591979688312, 879.3020219519976),
            yscale=:log10)
Plots.savefig(p3,"data/output/figures/final/projected-lambda1-linear.pdf")




#### COMBINING PLOTS
figs = Plots.plot(p3,p2,p1,
            layout=(1,3),size=(1400,400),
            bottom_margin=8Measures.mm,
            link=:y,
            yscale=:log10,
)
ylims(figs)
#visual spacing
Plots.plot!(figs[1],left_margin=12Measures.mm)
Plots.plot!(figs[2],left_margin=5Measures.mm)
Plots.plot!(figs[3],left_margin=5Measures.mm)
# Plots.plot!(figs[3],right_margin=8Measures.mm)
#link yaxes 
Plots.plot!(figs,link=:y)
Plots.plot!(figs,guidefontsize=16,tickfontsize=15,titlefontsize=20,
    bottom_margin=12Measures.mm)

# Plots.plot!(figs,background_color_inside=:gray98,top_margin=3Measures.mm)
#increase resolution
Plots.plot!(figs,dpi=300)

Plots.savefig(figs,"data/output/figures/final/projected-lambda1.pdf")
Plots.savefig(figs,"data/output/figures/poster/projected-lambda1.svg")