# for standalone running. loads in aggregated_data dictionaries  
# include("epidemic-figures-utils.jl")

# aggregated_data2 = load_aggregated_data("aggregated_data","aggregated-sirs-output-scratch-v2.json")

b = 1e-1
g = 5e-2 
hyper_beta = "linear"

d_vals = [1e-2,2e-2,3e-2]
for (ind, (d, tmp_data)) in enumerate(zip(d_vals, (aggregated_data2, aggregated_data3, aggregated_data3)))
    tmp = new_infections_data(tmp_data, "-50000-2", b, g, d, hyper_beta)
    qinfo = mapslices(x -> quantile(x, [0, 0.5, 1]), tmp, dims=1)'
        
    lower = qinfo[:, 1]
    middle = qinfo[:, 2]
    upper = qinfo[:, 3]

    fig = Plots.plot(ALPHA_VALS, middle, 
         ribbon = (middle - lower, upper - middle),
         fillalpha = 0.3,
         linewidth = 2,
         color = :blue,
         xlabel = L"\alpha",
         ylabel = "Average Total Infections",
         title = L"\delta=%$d",
         legend = false,
         size = (400, 300)
    )
    display(fig)
    Plots.savefig(fig,"data/output/figures/final/hyperedge-effects-p2-$ind.pdf")               
end

