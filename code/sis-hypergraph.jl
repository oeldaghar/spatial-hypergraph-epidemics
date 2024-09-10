#hypergraph represented as a sequence of hyperedges 




using DataStructures
using Random
using SparseArrays
# data structure 
mutable struct HyperSISData
    beta::Float64 #uniform infection probability
    gamma::Float64 #uniform recovery probability
    log1beta::Float64
    log1gamma::Float64
    hyperNeighbors::Vector{Vector{Vector{Int}}} #adjacent hyperneighbors
    itime::Vector{Int} # infection time
    rtime::Vector{Int} # last recovery time (I->S)
    snodes::Vector{Bool} # still suceptible nodes
    etimes::PriorityQueue #event times 
end
"""
    _clear!(E::SISData)

TBW
"""
function _clear!(E::SISData)
    fill!(E.itime, typemax(Int))
    fill!(E.rtime, typemax(Int))
    fill!(E.snodes, true)
    empty!(E.etimes)
end
# initialization
"""
    get_hyperneighbors(hedges::Vector{Vector{Int}})

TBW
"""
function get_hyperneighbors(hedges::Vector{Vector{Int}})
    #preprocess hedges 
    sort!.(hedges)
    sort!(hedges)
    unique!(hedges)
    #build ego-centric hyper_neighbors
    nnodes = maximum(maximum.(hedges))
    hyper_neighbors = [Vector{Vector{Int}}() for i=1:nnodes]
    for edge in hedges
        for (i,node) in enumerate(edge)
            push!(hyper_neighbors[node],vcat(edge[1:i-1],edge[i+1:end]))        
        end
    end
    return hyper_neighbors
end
function HyperSISData(hedges::Vector{Vector{Int}},
                        beta::Float64,
                        gamma::Float64)
    nnodes = size(A,1)
    @assert(0<beta<1, "must have 0<beta<1")

    itime = Vector{Int}(undef, nnodes)
    rtime = Vector{Int}(undef, nnodes)
    snodes = Vector{Bool}(undef, nnodes)
    etimes = PriorityQueue()
    hyper_neighbors = get_hyperneighbors(hedges)
    E = SISData(beta, gamma, 1/log1p(-beta), 1/log1p(-gamma),
            hyper_neighbors, itime, rtime, snodes, etimes)
    _clear!(E)
    return E
end
# infection event
#TODO figure out how to write infection event 
    #move node to i
    #for edge in hyper_neighbors
        #compute beta_h = f(|h|,|I⋂h|).
        # for susceptible node in edge 
            #attempt infection w/ beta_h 
    
# recovery event 
#TODO write this function
    #move node to S 
    #recompute infections from hyper neighbors

# main epidemic loop 
#TODO write this function

get_hyperneighbors(hedges)

