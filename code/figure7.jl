using JSON
using Plots
using Measures
using LaTeXStrings

function plot_data(trials)
    eigvals = reduce(hcat, trials) # now we have alpha as row index and trial as column index
    p = Plots.plot(xlabel="α", ylabel="Number of edges")  
    # show the full data as a ribbon
    minl = minimum.(eachrow(eigvals))
    maxl = maximum.(eachrow(eigvals))
    mid = (minl .+ maxl)/2
    upper = maxl .- mid
    lower = mid .- minl
    @show minl, maxl
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
projected_trials_squared = JSON.parsefile("data/output/spatial_graph_projected_lambda1_50000_2_squared.json")["projected_trials"]
projected_trials_sqrt = JSON.parsefile("data/output/spatial_graph_projected_lambda1_50000_2_sqrt.json")["projected_trials"]
projected_trials_linear = JSON.parsefile("data/output/spatial_graph_projected_lambda1_50000_2_linear.json")["projected_trials"]

projected_trials_squared = Vector{Vector{Float64}}(projected_trials_squared)
projected_trials_sqrt = Vector{Vector{Float64}}(projected_trials_sqrt)
projected_trials_linear = Vector{Vector{Float64}}(projected_trials_linear)

#standalone plots
p1 = plot_data(projected_trials_squared)
ticks = yticks(p1)[1]
Plots.plot!(p1, ylabel="",framestyle=:grid,
            yticks=(ticks[1],fill("",length(ticks[2]))))
Plots.plot!(p1,title="g(m)=m²",
    titlefont=font("Helvetica Bold", 14))

p2 = plot_data(projected_trials_sqrt)
ticks = yticks(p2)[1]
Plots.plot!(p2, ylabel="",framestyle=:grid,
            yticks=(ticks[1],fill("",length(ticks[2]))))
Plots.plot!(p2,title="g(m)=sqrt(m)",
    titlefont=font("Helvetica Bold", 14))

p3 = plot_data(projected_trials_linear)
Plots.plot!(p3, ylabel="Weighted Projected λ₁",framestyle=:grid,
            title="g(m)=m",
            titlefont=font("Helvetica Bold", 14))

#### COMBINING PLOTS
figs = Plots.plot(p3,p2,p1,
            layout=(1,3),size=(1400,400),
            bottom_margin=8Measures.mm,
)
#visual spacing
Plots.plot!(figs[1],left_margin=12Measures.mm)
Plots.plot!(figs[3],right_margin=8Measures.mm)
#link yaxes 
Plots.plot!(figs,link=:y)
Plots.plot!(figs,guidefontsize=18,tickfontsize=15,titlefontsize=20,
    bottom_margin=12Measures.mm)
Plots.plot!(figs,background_color_inside=:gray98,top_margin=3Measures.mm)
#increase resolution
Plots.plot!(figs,dpi=300)

Plots.savefig(figs,"data/output/figures/final/projected-lambda1.pdf")
Plots.savefig(figs,"data/output/figures/poster/projected-lambda1.svg")