#testing for defining experiment scope 
include("../code/sirs.jl")
include("../code/data-io.jl")

fnames = filter!(x->endswith(x,".txt"),readdir("data/hypergraphs/"))
fnames = get_fnames(get_gname(fnames[end]))

# beta,gamma,delta,exo,tmax = 1e-2,1e-2,1/365,5/1e6,365*10
beta,gamma,delta,exo,tmax = 0.9,1e-1,1/100,5/1e6,365*10

function single_epidemic_comparison(fnames;beta=3e-2,gamma=5e-2,delta=1/365,exo=5/10000,tmax=10*365)
    hedges = read_hedges(fnames[1])
    n = maximum(x->maximum(x),hedges)
    seed_node = rand(1:n)
    
    epidemic_results = []
    #epidemic simulations 
    for fname in fnames 
        hedges = read_hedges(fname)
        hyper_infections = first(sirs_hyper(hedges,seed_node,beta,gamma,delta,exo,tmax,hyper_beta_func="linear"))
        push!(epidemic_results,hyper_infections)
    end
    
    pairwise_infections = epidemic_results[1]
    hyper_infections = epidemic_results[end]

    
    computed_avgs = map(x->mean(vec(x)[end-365:end]),epidemic_results)
    print(typeof(computed_avgs))
    #figure layout 
    #put figures side by side 
    fig = Figure(resolution=(1000, 400))
    left_column = fig[1, 1] = GridLayout()
    right_column = fig[1, 2] = GridLayout()
    # Create axes in each column
    ax1 = Axis(left_column[1, 1], 
                title="HigherOrder Model",
                xlabel="Time",
                ylabel="Total Infections")
    ax2 = Axis(right_column[1, 1],
                title="Pairwise Model",
                xlabel="Time",
                ylabel="Total Infections")
    # Plot data in each axis
    lines!(ax1, hyper_infections, color=:blue)
    lines!(ax2, pairwise_infections, color=:red)
    linkaxes!(ax1, ax2)
    newfig = Figure(resolution=(800,2000))
    axes = []
    for i=1:lastindex(epidemic_results)
        alpha = _get_alpha(fnames[i])
        ax = Axis(newfig[i,1],title="alpha=$(alpha), time_lagged_avg=$(computed_avgs[i])")
        lines!(ax,epidemic_results[i])
        push!(axes,ax)
    end
    linkyaxes!(axes...)
    return fig,ax1,ax2,newfig
end

hedges = read_hedges(fnames[end])
#remove self loops 
sum(map(x->length(x)-length(unique(x)),hedges))

hedges[end]

tmp = _hyperedge_list_to_neighbor_list(hedges)
tmp[4992]



fig,ax1,ax2,newfig = single_epidemic_comparison(fnames,beta=beta,gamma=gamma,delta=delta,exo=exo,tmax=tmax)
fig
newfig
#run the above a few times. it seems that the epidemic peak is actually lower for the higherorder model despite what you would think (should be higher). perhaps analyze effects in terms of initial epidemic and then the cyclic/steady state analysis. I need to consider that this may be a bug..

xlims!(ax1,200,365*10)
xlims!(ax2,200,365*10)
ylims!(ax1,0,500)
ylims!(ax2,0,500)
display(fig)

#would have to dive deeper into the noise here. at a first glance, pairwise MIGHT be giving slightly higher steady state with more cyclic/oscillatory behavior but would have to scale up and quantify
newfig



fname = fnames[1]
hedges = read_hedges("$fname")


sum(map(x->length(x)-length(unique(x)),hedges))

699/length(hedges)

sirs_hyper(hedges,2)

Random.seed!(7)
d = 2
n = 10000
X = rand(d,n) #d = dimension, n = num nodes  
degreedist=LogNormal(log(3),1)
degs = rand(degreedist,n)
degs = min.(ceil.(Int,degs),n-1)

T = BallTree(X)
# form the edges for sparse
edges = Vector{Int}[]
# idxs, dists = knn(T, X[:,i], deg+1)
radfunc=(dist,deg) -> dist/sqrt(deg)
# Random.seed!(7)
for i=1:n #X
    deg = degs[i]
    idxs, dists = knn(T, X[:,i], deg+1)
    if deg > 1 
        maxdist = maximum(dists) 
        pts = @view X[:,idxs]
        rad = radfunc(maxdist,deg)
        # if !(rad ≈ maxdist)
        #   println("rad: ", rad, " maxdist: ", maxdist)
        # end
        clusters = dbscan(pts, rad)
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
    else
        # only one vertex! 
        push!(edges, [i,idxs[2]])
    end
end

i1 = findfirst(map(x->length(x)>length(unique(x)),edges))

edges[i1]

deg = degs[edges[i1][1]]
node = edges[i1][1]
idxs, dists = knn(T, X[:,node], 10)