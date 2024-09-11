#testing for defining experiment scope 
include("../code/sirs.jl")

fnames = filter!(x->endswith(x,".txt"),readdir("data/hypergraphs/"))
fname = fnames[1]
hedges = read_hedges("data/hypergraphs/$fname")
pairwise_edges = project(hedges)

beta,gamma,delta,exo,tmax = 3e-2,5e-2,1/365,15/1e6,10*365

function single_epidemic_comparison(hedges;beta=3e-2,gamma=5e-2,delta=1/365,exo=5/10000,tmax=10*365)
    pairwise_edges = project(hedges)
    n = maximum(x->maximum(x),hedges)
    seed_node = rand(1:n)
    #epidemic simulations 
    hyper_infections = sirs_hyper(hedges,seed_node,beta,gamma,delta,exo,tmax)
    pairwise_infections = sirs(pairwise_edges, seed_node, beta,gamma,delta,exo,tmax)

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
    # Display the figure
    return fig,ax1,ax2 
end

fig,ax1,ax2 = single_epidemic_comparison(hedges,beta=beta,gamma=gamma,delta=delta,exo=exo,tmax=tmax)
display(fig)
#run the above a few times. it seems that the epidemic peak is actually lower for the higherorder model despite what you would think (should be higher). perhaps analyze effects in terms of initial epidemic and then the cyclic/steady state analysis. I need to consider that this may be a bug..

xlims!(ax1,200,365*10)
xlims!(ax2,200,365*10)
ylims!(ax1,0,500)
ylims!(ax2,0,500)
display(fig)

#would have to dive deeper into the noise here. at a first glance, pairwise MIGHT be giving slightly higher steady state with more cyclic/oscillatory behavior but would have to scale up and quantify