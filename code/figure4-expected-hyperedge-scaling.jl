using Distributed
# Distributed.addprocs(10)

@everywhere using Distributions, Clustering, NearestNeighbors, Random, Plots, LaTeXStrings, Colors
@everywhere using StatsBase
@everywhere using ProgressMeter
@everywhere using MatrixNetworks, Graphs, SparseArrays
using Measures

@everywhere function hypergraph_edges(X,deg_list;radfunc=(dist,deg) -> dist/sqrt(deg))
    T = BallTree(X)
    edges = Vector{Int}[]
    dim = size(X,1)
    for i in axes(X,2)
        # deg = min(ceil(Int,rand(degreedist)),n-1)
        deg = deg_list[i]
        idxs, dists = knn(T, X[:,i], deg+1)
        if deg > 1 
            maxdist = maximum(dists) 
            pts = @view X[:,idxs]
            rad = radfunc(maxdist,deg,dim)
            clusters = dbscan(pts, rad).clusters
            for c in clusters
                e = [i]
                for v in c.core_indices
                    if idxs[v] != i
                        push!(e, idxs[v])
                    end
                end
                for v in c.boundary_indices
                    if idxs[v] != i 
                        push!(e, idxs[v])
                    end
                end
                if length(e) > 1
                    push!(edges, e)
                end
            end
        # else
            # # only one vertex! 
            # push!(edges, [i,idxs[2]])
        end 
    end 
    return edges
end

# make sure to update plot labels based on this function
@everywhere function get_func(alpha) 
    if alpha <= 1
        return (d,deg,dim) -> alpha^(1/dim)*d/(deg)^(1/(2*dim))
        # return (d,deg,dim) -> alpha*d/(deg)^(1/(2))
    else
        return (d,deg,dim) -> d/(deg)^(1/(2*dim)) + (alpha-1)*(d-d/(deg)^(1/(2*dim)))
        # return (d,deg,dim) -> d/(deg)^(1/(2)) + (alpha-1)*(d-d/(deg)^(1/(2)))
    end
end 

Random.seed!(1234)

# @everywhere d = 2
# @everywhere n = 1000
# @everywhere alphas = range(0,1,length=10) # first portion 

@everywhere function run_trial(n,d,alphas)
    X = rand(d,n)

    deg_list = zeros(Int, n)
    degreedist = LogNormal(log(3),1)
    for i = 1:n
        deg_list[i] = min(ceil(Int,rand(degreedist)),n-1)
    end

    exp_sqrt_deg = mean(sqrt.(deg_list))

    graphs = pmap(alpha -> hypergraph_edges(X, deg_list;radfunc=get_func(alpha)), alphas)
    avg_hsize = pmap(x->mean(length.(x)), graphs)
    # alpha^d * E(\sqrt(k_v))
    rhs = exp_sqrt_deg
    # depends on the function for epsilon
    lowerbounds = 1.0 .+ rhs.*(alphas)

    return avg_hsize, lowerbounds
end
# trials = @showprogress map(x->run_trial(n,d,alphas), 1:20)

function make_plot(n,d,alphas,ntrials)
    trials = @showprogress map(x->run_trial(n,d,alphas), 1:ntrials)
    hsizes = reduce(hcat, first.(trials))
    lowerbounds = reduce(hcat, map(x->x[2],trials))
    # compute averages 
    hsizes = vec(mean(hsizes,dims=2))
    lowerbounds = vec(mean(lowerbounds,dims=2))

    plt = Plots.plot()
    Plots.plot!(plt, size=(500,500))
    Plots.plot!(plt, xlabel=L"\alpha")
    Plots.plot!(plt, alphas, hsizes, label="", color=:blue)
    Plots.plot!(plt, alphas, lowerbounds, label="", color=:grey)
    Plots.plot!(plt,plot_title="n=$n, d=$d, \neps = alpha^1/d r_v/k_v^(1/2d)")
    Plots.plot!(plt,annotation=(0.75,2.15,L"\alpha E(\sqrt{k_v})"))
    Plots.plot!(plt,annotation=(0.25,2.75,L"E( | h |)"),c=:blue)
    display(plt)
    return plt 
end


n = 2500
d = 25
alphas = range(0,1,length=10)


make_plot(5000,2,alphas,20)

make_plot(5000,3,alphas,20)

make_plot(5000,5,alphas,20)

make_plot(5000,25,alphas,10)

