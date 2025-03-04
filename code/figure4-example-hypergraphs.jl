### plotting graphs 
Random.seed!(11235)
n = 250
d = 2
X = rand(d,n)
deg_list = zeros(Int, n)
degreedist = LogNormal(log(3),1)
for i = 1:n
    deg_list[i] = min(ceil(Int,rand(degreedist)),n-1)
end

# make graphs and convert to biadjacency matrices 
graphs = pmap(alpha -> hypergraph_edges(X, deg_list;radfunc=get_func(alpha)), range(0,2,3))

function plot_hypergraph(X,hedges,nodes; a=0.1)
    plt = Plots.plot(leg=false, size=(500,500))
    for edge in hedges
        if length(edge) > 2 
            local b = [Ball2(X[:, v], 0.0005) for v in edge]
            local c = ConvexHullArray(b)
            Plots.plot!(plt, c, 1e-3, alpha=a, c=:blue)
        else
            Plots.plot!(plt, X[1,edge], X[2,edge], c=:blue,alpha=0.35)
        end
    end
    # xs = [X[1,v] for v in nodes]
    # ys = [X[2,v] for v in nodes]
    # Plots.scatter!(plt,xs,ys,leg=false,markerstrokewidth=0,markersize=12,markercolor=:blue)
    return plt
end

figs = []
for (ind,h) in enumerate(graphs)
    f = plot_hypergraph(X,h,1:n,a=0.1)
    Plots.plot!(f,framestyle=:none,
                margins=-20Plots.mm,
                )
    Plots.scatter!(f,X[1,:],X[2,:],leg=false,
        markerstrokewidth=0,
        markercolor=:black,
        markersize=3,
    )
    display(f)    
    push!(figs,f)
    Plots.savefig(f,"data/output/figures/final/graph-demo-$ind.pdf")
end