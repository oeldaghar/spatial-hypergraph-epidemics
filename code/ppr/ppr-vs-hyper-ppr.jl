using Distributed
# Distributed.addprocs(10)
@everywhere include("../spatial-hypergraph.jl")
@everywhere using Distributions, Random
@everywhere using StatsBase
@everywhere using ProgressMeter
@everywhere using MatrixNetworks, Graphs, SparseArrays
@everywhere using SparseArrays
using LazySets
# Plotting 
using Plots, LaTeXStrings, Colors
using Measures

# uses Meng's code from https://github.com/MengLiuPurdue/LHQD/tree/main
include("LHQD/local-hyper.jl")
include("LHQD/common.jl")
### EXAMPLE FROM MENG DOCS
# H should be biadjacency matrix of the form E X V (SparseCSC)
# G = LH.graph(H,1.0) # a hypergraph object with delta=1.0 in its cut function
# q = 2.0
# L = LH.loss_type(q) # the loss-type, this is a 2-norm)
# kappa = 0.01 # value of kappa (sparsity regularization)
# gamma = 0.1 # value of gamma (regularization on seed) 
# rho = 0.5 # value of rho (KKT apprx-val)
# S = [1] # seed set
# x,r,iter = LH.lh_diffusion(G,S,gamma,kappa,rho,L,max_iters=10000)
# cond,cluster = hyper_sweepcut(G.H,x,G.deg,G.delta,0.0,G.order)
### END EXAMPLE

# hedges to biadjacency E X V as SparseCSC
function hedges_to_biadj(hedges,n)
    # hedges[hyperedge_id] = [nodes in hyperedge] 
    nedges = lastindex(hedges)
    
    nodes = Vector{Int}()
    edges = Vector{Int}()
    for (row,h) in enumerate(hedges)
        append!(nodes,h)
        append!(edges,row*ones(lastindex(h)))
    end
    return sparse(edges,nodes,ones(lastindex(nodes)),lastindex(hedges),n)
end

function pairwise_proj(hedges)
    ei = Vector{Int}()
    ej = Vector{Int}()
    for h in hedges
        for node1 in h 
            for node2 in h 
                push!(ei,node1)
                push!(ej,node2)
            end
        end
    end
    nnodes = maximum([maximum(ei),maximum(ej)])
    A = sparse(ei,ej,ones(lastindex(ei)),nnodes,nnodes)
    for i=1:nnodes
        A[i,i]=0
    end
    dropzeros!(A)
    return A
end

# take hypergraph 
# pairwise_proj 
# run ppr diffusions from same seed set 

n,d = 500,2
alphas = range(0,2,3)

# hypergraph ppr parameters 
q = 2.0
L = LH.loss_type(q) # the loss-type, this is a 2-norm)
kappa = 0.0001 # value of kappa (sparsity regularization)
gamma = 1.0 # value of gamma (regularization on seed) 
rho = 0.5 # value of rho (KKT apprx-val)
# S = [rand(1:n)] # seed set

Random.seed!(11)
X = rand(d,n)

_,S = findmin(vec(mapslices(xy->(xy[1]-1)^2+(xy[2]-1)^2,X,dims=1)))
S = [S]


deg_list = zeros(Int, n)
degreedist = LogNormal(log(3),1)
for i = 1:n
    deg_list[i] = min(ceil(Int,rand(degreedist)),n-1)
end

# make graphs and convert to biadjacency matrices 
graphs = pmap(alpha -> hypergraph_edges(X, deg_list;radfunc=get_func(alpha)), alphas)

# hypergraph ppr
function hypergraph_ppr(H,S)
    # H = hedges_to_biadj(hedges,n) # biadjacency rep
    G = LH.graph(H,1.0) # a hypergraph object with delta=1.0 in its cut function
    x,r,iter = LH.lh_diffusion(G,S,gamma,kappa,rho,L,max_iters=10000)
    # cond,cluster = hyper_sweepcut(G.H,x,G.deg,G.delta,0.0,G.order)
    return x 
end

# hedges = graphs[end]
# H = hedges_to_biadj(hedges,n)
# ppr_soln = hypergraph_ppr(H,S)
# newx = log10.(1e-16 .+ ppr_soln)

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
    # xs = [X[1,v] for v in nodes]
    # ys = [X[2,v] for v in nodes]
    # Plots.scatter!(plt,xs,ys,leg=false,markerstrokewidth=0,markersize=12,markercolor=:blue)
    return plt
end

# plt = plot_hypergraph(X,hedges,1:n)
# Plots.scatter!(X[1,:],X[2,:],marker_z=newx)

function make_ppr_plot(hedges)
    H = hedges_to_biadj(hedges,n)
    ppr_soln = hypergraph_ppr(H,S)
    new_inds = findall(ppr_soln.>0)

    newx = log10.(1e-16 .+ ppr_soln)
    markersizes = 3 .*ones(Float64,lastindex(newx))
    markersizes[ppr_soln.>0] .= 8

    plt = plot_hypergraph(X,hedges,1:n)
    Plots.scatter!(plt,X[1,:],X[2,:],marker_z=newx,
            markerstrokewidth=0.5,
            markersize=markersizes,
            xlabel = "Spatial Coordinate 1",
            ylabel = "Spatial Coordinate 2",
            framestyle=:box,
            # xlims=(0.5,1.05),
            # ylims = (0.5,1.05),
            )
    Plots.scatter!(plt,X[1,new_inds],X[2,new_inds],marker_z=newx[new_inds],
        markerstrokewidth=0.5,
        markersize=markersizes[new_inds],
        xlabel = "Spatial Coordinate 1",
        ylabel = "Spatial Coordinate 2",
        framestyle=:box,
        # xlims=(0.5,1.05),
        # ylims = (0.5,1.05),
        )
    return plt 
end

figs = map(x->make_ppr_plot(x),graphs)
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



ppr_soln = map(x->hypergraph_ppr(hedges_to_biadj(x,n),S),graphs)
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