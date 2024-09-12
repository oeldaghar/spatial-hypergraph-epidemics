#=
    this code was modified from code sent by david. 
    not sure why the discrete event simulation is slower than this..

    update: it looks like there's a ton of overhead in the PriorityQueue from DataStructures.. will figure out later
=#

using Random, MatrixNetworks
include("hypergraph-tools.jl")
include("spatial-hypergraph.jl")

using ProgressMeter
using Distributed

function _edge_list_to_neighbor_list(edges)
    n = maximum(x->max(x[1],x[2]),edges) # find the number of nodes
    neighbors = [Int[] for i in 1:n] # initialize the neighbor list
    for edge in edges
        push!(neighbors[edge[1]], edge[2]) # add the second node to the list of neighbors of the first node
        push!(neighbors[edge[2]], edge[1]) # add the first node to the list of neighbors of the second node
    end
    return neighbors
end
# create an enum for the possible states
@enum NodeState StateSucceptible=1 StateInfected=2 StateRecovered=3
function _update_state!(neighborlist, currstate, newstate, beta, gamma, delta, exo=5/100000)
    copy!(newstate, currstate) # copy the current state to the new state
    for i in eachindex(neighborlist)
        if rand() < exo # exogenous infection #5 per 100k cases per unit time
            newstate[i] = StateInfected
            continue 
        end
        if currstate[i] == StateSucceptible
            #check neighbors 
            for j in neighborlist[i]
                if currstate[j] == StateInfected
                    if rand() < beta
                        newstate[i] = StateInfected
                        break # no need to check the other neighbors if one is infected
                    end
                end
            end
        elseif currstate[i] == StateInfected
            if rand() < gamma
                newstate[i] = StateRecovered
            end
        elseif currstate[i] == StateRecovered
            if rand() < delta
                newstate[i] = StateSucceptible
            end
        end
    end
end
function sirs(edges::Vector{T}, seednodes::Vector{Int}, beta=3e-2, gamma=5e-2, delta=5e-2, exo=5/100000, tmax=365*10) where T
    neighlist = _edge_list_to_neighbor_list(edges)
    n = length(neighlist)
    
    state = Vector{NodeState}(undef, n)
    fill!(state, StateSucceptible)
    for node in seednodes
        state[node] = StateInfected
    end
    newstate = Vector{NodeState}(undef, n)

    ninfected = zeros(Int, tmax)

    @showprogress for t in 1:tmax
        ninfected[t] = count(x->x==StateInfected, state)
        _update_state!(neighlist, state, newstate, beta, gamma, delta, exo)
        state, newstate = newstate, state # swap... 
    end
    return ninfected
end 
sirs(edges, seednode::Int, beta=3e-2, gamma=5e-2, delta=5e-2, exo=5/100000, tmax=365*10) = sirs(edges, [seednode], beta, gamma, delta, exo, tmax)

#hypergraph version 
"""
    _hyperedge_list_to_neighbor_list(hedges)

convert hyper edges in to neighbor list. if [i,other nodes] is a hyperedge,
then hyper_neighbors[i] contains [other nodes].
"""
function _hyperedge_list_to_neighbor_list(hedges)
    #preprocess hedges 
    sort!.(hedges)
    sort!(hedges)
    unique!(hedges)
    
    n = maximum(x->maximum(x),hedges) 
    
    #build ego-centric hyper_neighbors with weights.
    # TODO add weights so that 
    # neighs[i] = [(h,w_h): h is ego hyperedge for i and w_h is weight]
    hyper_neighbors = [Vector{Vector{Int}}() for i=1:n]
    for edge in hedges
        for (i,node) in enumerate(edge)
            push!(hyper_neighbors[node],vcat(edge[1:i-1],edge[i+1:end]))        
        end
    end
    return hyper_neighbors
end
#TODO cache info on hyperedge stats (number of infected nodes)
#try caching at beginning of iteration and do a lookup
function _update_state_hyper!(hyper_neighborlist, currstate, newstate, beta, gamma, delta, exo=5/100000)
    copy!(newstate, currstate) # copy the current state to the new state
    for i in eachindex(hyper_neighborlist)
        if rand() < exo # exogenous infection #5 per 100k cases per unit time
            newstate[i] = StateInfected
            continue 
        end
        if currstate[i] == StateSucceptible
            #check hyper neighbors 
            for i_ego_edge in hyper_neighborlist[i] #hedge = [i,other nodes], then ego_edge = [other nodes]
                #measure fraction of nodes infected 
                infected_ego_size = 0
                for j in i_ego_edge
                    if currstate[j] == StateInfected
                        infected_ego_size+=1
                    end
                end
                #define hyperedge contribution
                # hyper_beta = f(|h|,|h ⋂ I|)... figure out later
                if infected_ego_size>0
                    # hyper_beta = beta #change later..
                    hyper_beta = 1-(1-beta)^infected_ego_size 
                    hyper_beta/=lastindex(i_ego_edge)
                    # hyper_beta/=lastindex(i_ego_edge)
                    #alter this probability downwards a term like sqrt(|h|). something about air to room ratio and airflow. check details.
                    if rand() < hyper_beta
                        newstate[i] = StateInfected
                        break # no need to check the other neighbors if one is infected
                    end
                end
            end
        elseif currstate[i] == StateInfected
            if rand() < gamma
                newstate[i] = StateRecovered
            end
        elseif currstate[i] == StateRecovered
            if rand() < delta
                newstate[i] = StateSucceptible
            end
        end
    end
end
function sirs_hyper(hedges, seednodes::Vector{Int}, beta=1e-2, gamma=5e-2, delta=1/365, exo=5/100000, tmax=365*10)
    neighlist = _hyperedge_list_to_neighbor_list(hedges)
    n = length(neighlist)
    
    state = Vector{NodeState}(undef, n)
    fill!(state, StateSucceptible)
    for node in seednodes
        state[node] = StateInfected
    end
    newstate = Vector{NodeState}(undef, n)

    ninfected = zeros(Int, tmax)

    @showprogress for t in 1:tmax
        ninfected[t] = count(x->x==StateInfected, state)
        _update_state_hyper!(neighlist, state, newstate, beta, gamma, delta, exo)
        state, newstate = newstate, state # swap... 
    end
    return ninfected,state
end 
sirs_hyper(hedges, seednode::Int, beta=3e-2, gamma=5e-2, delta=5e-2, exo=5/100000, tmax=365*10) = sirs_hyper(hedges, [seednode], beta, gamma, delta, exo, tmax)

#TODO update with previous changes 
### parallel version of the code
function parallel_sirs(edges, seed_nodes;
                        beta=3e-2, gamma=5e-2, delta=1/365, exo=5/100000, tmax=365*10)
    println("USING $(lastindex(workers())) WORKERS for $(lastindex(seed_nodes)) NODES")
    #cache worker information
    @everywhere edges = $edges 
    #temporary.. if we want to vary params, wrap this in an iterator and modify the below helper function
    @everywhere beta,gamma,delta,exo,tmax = $beta,$gamma,$delta,$exo,$tmax
    @everywhere cached_edges_function(seednode) = first(sirs(edges,seednode,beta,gamma,delta,exo,tmax)) #get infection results but toss out states
    #main pmap function call
    infection_results = @showprogress pmap(x->cached_edges_function(x),seed_nodes)
    #clean up
    return infection_results
end 

function parallel_sirs_hyper(hedges, seed_nodes;
                        beta=3e-2, gamma=5e-2, delta=1/365, exo=5/100000, tmax=365*10)
    println("USING $(lastindex(workers())) WORKERS for $(lastindex(seed_nodes)) NODES")
    #cache worker information
    @everywhere hedges = $hedges 
    #temporary.. if we want to vary params, wrap this in an iterator and modify the below helper function
    @everywhere beta,gamma,delta,exo,tmax = $beta,$gamma,$delta,$exo,$tmax
    @everywhere cached_edges_function(seednode) = first(sirs_hyper(hedges,seednode,beta,gamma,delta,exo,tmax)) 
    #main pmap function call
    infection_results = @showprogress pmap(x->cached_edges_function(x),seed_nodes)
    #clean up
    return infection_results
end

