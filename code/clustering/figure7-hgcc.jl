using Plots, LaTeXStrings, Colors
using Measures
using JSON 

n = 10000
println("Reading data...")
loaded_data = open("data/output/hgcc_data-n_$n.json", "r") do file
    JSON.parse(file)
end

vrow1 = [
    loaded_data["($n, 2, \"global_cc_g_mat\")"],
    loaded_data["($n, 2, \"local_cc_g_mat\")"],
    loaded_data["($n, 2, \"global_cc_hg_mat\")"],
    loaded_data["($n, 2, \"local_cc_mg_mat\")"]
]

vrow2 = [
    loaded_data["($n, 5, \"global_cc_g_mat\")"],
    loaded_data["($n, 5, \"local_cc_g_mat\")"],
    loaded_data["($n, 5, \"global_cc_hg_mat\")"],
    loaded_data["($n, 5, \"local_cc_mg_mat\")"]
]

vrow3 = [
    loaded_data["($n, 10, \"global_cc_g_mat\")"],
    loaded_data["($n, 10, \"local_cc_g_mat\")"],
    loaded_data["($n, 10, \"global_cc_hg_mat\")"],
    loaded_data["($n, 10, \"local_cc_mg_mat\")"]
]

vrow1 = map(x->reduce(hcat,x),vrow1)
vrow2 = map(x->reduce(hcat,x),vrow2)
vrow3 = map(x->reduce(hcat,x),vrow3)

# testing 
# norm(norm.(vrow1.-row1))
# norm(norm.(vrow2.-row2))
# norm(norm.(vrow3.-row3))

# PLOT 1 - quantiles
function quantile_plot(data_mat,alphas)
    q_info = mapslices(x->quantile(x,[0.0,0.5,1.0]),data_mat,dims=2)
    q_lower = q_info[:,2].-q_info[:,1]
    q_upper = q_info[:,3].-q_info[:,2]
    f = Plots.plot(alphas, q_info[:,2], ribbon=(q_lower,q_upper), 
                leg=false, 
                # xlabel=L"\alpha", 
                # ylabel="Pairwise CC", 
                guidefontsize = 14,
                linewidth=2, 
                marker=:circle, 
                markerstrokewidth=0, 
                ylims=(0,1),
    )   
    return f 
end

function make_cc_fig(row_data)
    global_cc_g_mat, local_cc_g_mat ,global_cc_hg_mat, local_cc_mg_mat = row_data

    f1 = quantile_plot(global_cc_g_mat,alphas)
    Plots.plot!(f1,
                xlabel=L"\alpha", 
                ylabel="Global Pairwise CC", 
            )

    f2 = quantile_plot(local_cc_g_mat,alphas)
    Plots.plot!(f2,
                xlabel=L"\alpha", 
                ylabel="Local Pairwise CC", 
                ylims = (0.2,0.85)
            )

    f3 = quantile_plot(global_cc_hg_mat,alphas)
    Plots.plot!(f3,
                xlabel=L"\alpha", 
                ylabel=L"CC_4", 
                ylims=(0,0.35)
            )

    f4 = quantile_plot(local_cc_mg_mat,alphas)
    Plots.plot!(f4,
                xlabel=L"\alpha", 
                ylabel="Multigraph Local CC", 
                ylims=(0.15,0.5)
            )
    
    plt = Plots.plot(f1,f2,f4,f3,layout=(1,4), margin=6*Plots.mm,size=(1700,400))
                                        
    #touch up margins
    Plots.plot!(plt[2],left_margin=2Measures.mm,bottom_margin=5Measures.mm)
    Plots.plot!(plt[1],left_margin=10Measures.mm)
    return plt 
end


f1 = make_cc_fig(vrow1)
Plots.plot!(f1,plot_title="n=$n d=2",plot_titlefontsize=24)
f2 = make_cc_fig(vrow2)
Plots.plot!(f2,plot_title="n=$n d=5",plot_titlefontsize=24)
f3 = make_cc_fig(vrow3)
Plots.plot!(f3,plot_title="n=$n d=10",plot_titlefontsize=24)

# put them all together
plt = Plots.plot(f1,f2,f3,layout=(3,1),size=(1400,1200),
        top_margin=-5Plots.mm,dpi=500)

# Plots.savefig(plt,"data/output/figures/testing/hypergraph-cc.png")
Plots.savefig(plt,"data/output/figures/final/hypergraph-cc.pdf")
