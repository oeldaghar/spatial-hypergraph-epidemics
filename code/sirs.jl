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
# replacing enum with smaller memory footprint based on profiling results.
const StateSucceptible = Int8(1)
const StateInfected = Int8(2)
const StateRecovered = Int8(3)
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
    
    state = Vector{Int8}(undef, n)
    fill!(state, StateSucceptible)
    for node in seednodes
        state[node] = StateInfected
    end
    newstate = Vector{Int8}(undef, n)

    ninfected = zeros(Int, tmax)

    @showprogress for t in 1:tmax
        ninfected[t] = count(x->x==StateInfected, state)
        _update_state!(neighlist, state, newstate, beta, gamma, delta, exo)
        state, newstate = newstate, state # swap... 
    end
    return ninfected,state
end 
sirs(edges, seednode::Int, beta=3e-2, gamma=5e-2, delta=5e-2, exo=5/100000, tmax=365*10) = sirs(edges, [seednode], beta, gamma, delta, exo, tmax)

# hypergraph version 
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

# struct to store ego_edge (hyper neighbors) and the weight (number infected nodes in hyperedge)
mutable struct Neighbor
    ego_edge::Vector{Int}
    weight::Int
end

# hyper_beta_func options: "linear", "sqrt", "squared", or "pairwise"

function hyperedge_normalization(hyper_beta_func::String, size::Int)
    if hyper_beta_func == "linear"
        return size
    elseif hyper_beta_func == "sqrt"
        return sqrt(size)
    elseif hyper_beta_func == "squared"
        return size^2
    elseif hyper_beta_func == "pairwise"
        return 1.0
    else
        throw(ArgumentError("hyper_beta_func should be one of 'linear', 'sqrt', 'squared', or 'pairwise'"))
    end
end
function _update_state_hyper!(
                                weighted_neighlist::Vector{Vector{Neighbor}},
                                currstate::Vector{Int8}, 
                                newstate::Vector{Int8}, 
                                beta::Float64, 
                                gamma::Float64, 
                                delta::Float64, 
                                exo::Float64=5/100000,
                                hyper_beta_func::String="linear",
                                ninfected_hyperedge::Dict{Int,Int}=Dict{Int,Int}(),
                                ntransmissions_hyperedge::Dict{Int,Int}=Dict{Int,Int}()
                            )
    copy!(newstate, currstate) # copy the current state to the new state)
    # caching this for slight optimizatioin
    beta_complement = 1-beta 
    @inbounds for i in eachindex(weighted_neighlist)
        if rand() < exo # exogenous infections
            newstate[i] = StateInfected
            continue 
        end
        if currstate[i] == StateSucceptible
            #check hyper neighbors 
            for hyper_neighbor in weighted_neighlist[i] #hedge = [i,other nodes], then ego_edge = [other nodes]
                i_ego_edge,infected_ego_size = hyper_neighbor.ego_edge,hyper_neighbor.weight
                #define hyperedge contribution
                # hyper_beta = f(|h|,|h ⋂ I|) .. depends on size of hyperedge and number of infected nodes 
                if infected_ego_size>0
                    # 1-(1-beta)^(num_infected_nodes) / hyperedge_normalization
                    hyper_beta = 1-beta_complement^infected_ego_size 
                    hyper_beta /= hyperedge_normalization(hyper_beta_func, lastindex(i_ego_edge))
                    if rand() < hyper_beta
                        newstate[i] = StateInfected

                        # record size of hyperedge that passed infection 
                        h_size = lastindex(i_ego_edge)+1
                        old_val = get(ntransmissions_hyperedge,h_size,0)
                        ntransmissions_hyperedge[h_size] = old_val+1
                        # break # no need to check the other neighbors if one is infected
                        # if we want to allow ties, we need to remove this 
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
    #update weights in weighted_hyper_neighlist 
    @inbounds for node=1:lastindex(weighted_neighlist)
        for hyper_neighbor in weighted_neighlist[node]
            ego_edge = hyper_neighbor.ego_edge
            s = 0 
            for i in ego_edge
                if newstate[i]==StateInfected
                    s+=1
                end
            end
            #update weight
            hyper_neighbor.weight = s

            # record number of infected nodes in binned infections 
            h_size = lastindex(ego_edge)+1
            old_val = get(ninfected_hyperedge,h_size,0)
            ninfected_hyperedge[h_size] = old_val+s+newstate[node] # sum_v (inf neghbors of v + v==infected)
        end
    end
end

function sirs_hyper(hedges::Vector{Vector{Int}}, seednodes::Vector{Int}, beta::Float64=1e-2, gamma::Float64=5e-2, delta::Float64=1/365, exo::Float64=5/100000, tmax::Int=365*10,hyper_beta_func::String="linear")
    # build adjacency list for hyper-neighbors 
    neighlist = _hyperedge_list_to_neighbor_list(hedges)
    # weighted neighbor list for number of infected nodes in each hyperedge 
    weighted_neighlist = Vector{Vector{Neighbor}}(undef, length(neighlist))
    for i in eachindex(neighlist)
        weighted_neighlist[i] = [Neighbor(y, 0) for y in neighlist[i]]
    end

    n = length(neighlist)
    # initialize infected nodes
    state = Vector{Int8}(undef, n)
    fill!(state, StateSucceptible)
    for node in seednodes
        state[node] = StateInfected
    end
    #update weights for infected nodes in each hyperedge 
    for node=1:lastindex(weighted_neighlist)
        for hyper_neighbor in weighted_neighlist[node]
            ego_edge = hyper_neighbor.ego_edge
            s = 0 
            for i in ego_edge
                if state[i]==StateInfected
                    s+=1
                end
            end
            #update weight
            hyper_neighbor.weight = s
        end
    end

    # initialize info 
    newstate = Vector{Int8}(undef, n)
    ninfected = zeros(Int, tmax)
    ninfected_hyperedge_size = Vector{Dict{Int,Int}}(undef,tmax) # number of infected nodes in a binned hyperedge 
    nhyperedge_transmissions = Vector{Dict{Int,Int}}(undef,tmax) # number of times a binned hyperedge passes an infection 
    @showprogress for t in 1:tmax # main loop 
        ninfected[t] = count(x->x==StateInfected, state)
        # intialize binned hyperedge stats 
        t_ninfected_hyperedge_size = Dict{Int,Int}()
        t_nhyperedge_transmissions = Dict{Int,Int}()
        # inplace updates 
        _update_state_hyper!(weighted_neighlist, state, newstate, beta, gamma, delta, exo, hyper_beta_func,
                                t_ninfected_hyperedge_size,
                                t_nhyperedge_transmissions)
        # update states 
        state, newstate = newstate, state # swap... 
        # cache binned hyperedge stats 
        ninfected_hyperedge_size[t] = Dict(key=>val for (key,val) in pairs(t_ninfected_hyperedge_size))
        nhyperedge_transmissions[t] = Dict(key=>val for (key,val) in pairs(t_nhyperedge_transmissions))
    end
    return ninfected,state,ninfected_hyperedge_size,nhyperedge_transmissions
end 
sirs_hyper(hedges::Vector{Vector{Int}}, seednode::Int, beta=3e-2, gamma=5e-2, delta=5e-2, exo=5/100000, tmax=365*10, hyper_beta_func::String="linear") = sirs_hyper(hedges, [seednode], beta, gamma, delta, exo, tmax,hyper_beta_func)

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
    infection_results = @showprogress pmap(x->compute_graph_stats(x),seed_nodes)
    #clean up
    return infection_results
end 
# function parallel_sirs_hyper(hedges, seed_nodes;
#                         beta=3e-2, gamma=5e-2, delta=1/365, exo=5/100000, tmax=365*10, hyper_beta_func="linear")
#     println("USING $(lastindex(workers())) WORKERS for $(lastindex(seed_nodes)) NODES")
#     #cache worker information
#     @everywhere hedges = $hedges 
#     #temporary.. if we want to vary params, wrap this in an iterator and modify the below helper function
#     @everywhere beta,gamma,delta,exo,tmax,hyper_beta_func = $beta,$gamma,$delta,$exo,$tmax,$hyper_beta_func
#     @everywhere cached_edges_function(seednode) = first(sirs_hyper(hedges,seednode,beta,gamma,delta,exo,tmax,hyper_beta_func)) 
#     #main pmap function call
#     infection_results = @showprogress pmap(x->cached_edges_function(x),seed_nodes)
#     #clean up
#     return infection_results
# end

# want to do all parameters for a single network. for epidemic parameters, can accept either Float64 or Vector{Float64}
# for seed_nodes, could be a vector or vector of vectors. 
function parallel_sirs_hyper(
    hedges,
    seed_nodes::Vector{Vector{S}} = [[1]];
    beta::Vector{T} = [5e-2], 
    gamma::Vector{T} = [5e-2], 
    delta::Vector{T} = [1/365], 
    exo::Vector{T} = [5/100000], 
    tmax::Vector{S} = [365*10], 
    hyper_beta_func::Vector{String} = ["linear"]
) where {T, S}

    #cache worker information
    @everywhere hedges = $hedges 
    # wrap other args into an iterator 
    epi_params = Base.Iterators.product(seed_nodes,beta,gamma,delta,exo,tmax,hyper_beta_func)
    # helper function for passing parameters. 
    @everywhere compute_graph_stats(seeds,beta,gamma,delta,exo,tmax,hyper_beta_func) = sirs_hyper(hedges,seeds,beta,gamma,delta,exo,tmax,hyper_beta_func)[[1,3,4]]
    
    println("USING $(lastindex(workers())) WORKERS for $(length(epi_params)) PARAMETERS")
    #main parallel function call
    infection_results = @showprogress pmap(x->compute_graph_stats(x...), epi_params)
    #clean up
    return infection_results, epi_params
end
