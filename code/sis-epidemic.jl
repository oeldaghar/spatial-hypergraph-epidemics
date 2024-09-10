module EventEpidemic
using DataStructures
using Random 
using SparseArrays
using ProgressMeter
mutable struct SISData
    beta::Float64 #uniform infection probability
    gamma::Float64 #uniform recovery probability
    log1beta::Float64
    log1gamma::Float64
    inNeighborsA::Vector{Vector{Int}} #adjlist representation of in-neighbors of A: (inneighbor, weight)
    outNeighborsA::Vector{Vector{Int}} #adjlist representation of out-neighbors of A: (outneighbor, weight)
    itime::Vector{Int} # infection time
    rtime::Vector{Int} # last recovery time (I->S)
    snodes::Vector{Bool} # still suceptible nodes
    etimes::PriorityQueue{Int,Int} #event times 
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
"""
    getneighbors(A::SparseMatrixCSC)

TBW
"""
function getneighbors(A::SparseMatrixCSC)
    inNeighbors = Vector{Vector{Int}}(undef,max(size(A)...))
    rowids = rowvals(A)
    @inbounds for i = 1:max(size(A)...)
        inNeighbors[i] = rowids[nzrange(A,i)]
    end
    return inNeighbors
end
"""
    SISData(A,beta::Float64,gamma::Float64)

TBW
"""
function SISData(A::SparseMatrixCSC,beta::Float64,gamma::Float64)
    nnodes = size(A,1)
    @assert(0<beta<1, "must have 0<beta<1")

    itime = Vector{Int}(undef, nnodes)
    rtime = Vector{Int}(undef, nnodes)
    snodes = Vector{Bool}(undef, nnodes)
    etimes = PriorityQueue{Int,Int}()
    inNeighborsA = getneighbors(A)
    outNeighborsA = getneighbors(sparse(A'))
    E = SISData(beta, gamma, 1/log1p(-beta), 1/log1p(-gamma),
            inNeighborsA, outNeighborsA, itime, rtime, snodes, etimes)
    _clear!(E)
    return E
end
"""
    _infect_node!(E::SISData, t::Int, n::Int)

TBW
"""
function _infect_node!(E::SISData, t::Int, n::Int)
    if E.snodes[n]
        E.snodes[n] = false
        etime = t # exposure time is now for node n
        ctime = etime+1 # forces nodes to be in E for 1 day
        rtime = ctime+(-randexp()*E.log1gamma) #rand(Geometric(E.gamma))
        rtime = ceil(Int,rtime)
        E.itime[n] = ctime
        E.rtime[n] = rtime 
        #schedule event for recovery  
        E.etimes[n] = rtime
        # handle neighbor updates 
        @inbounds for (j,v) in enumerate(E.outNeighborsA[n])
            if E.snodes[v] # not yet handled
                #dt = rand(Geometric(E.beta)) # number of interactions needed for an infection
                # rand(rng::AbstractRNG, d::Geometric) = floor(Int,-randexp(rng) / log1p(-d.p))
                dt = -randexp() * E.log1beta
                dti = floor(Int,dt)
                itimev = ctime+dti
                olditime = E.itime[v]
                if olditime<t || olditime==typemax(Int) #ignore old infection time
                    if itimev < rtime # then we may infect
                        E.itime[v] = itimev
                        E.etimes[v] = itimev # set the next event time for this node.
                    end
                end
            end 
        end #end neighbor loop
    end
end
"""
    _recovery(E,n)

TBW
"""
function _recovery!(E,t,n)
    #recover node
    if !(E.snodes[n])
        E.snodes[n] = true
        #reset itime 
        E.itime[n] = typemax(Int)
    end
    #fetch infectious neighbors and attempt to reinfect.
    @inbounds for v in E.inNeighborsA[n]
        if !(E.snodes[v])
            dt = -randexp() * E.log1beta
            dti = floor(Int,dt)
            curr_itime_n = t+dti
            if curr_itime_n < min(E.rtime[v],E.itime[n]) # then we may infect
                E.itime[n] = curr_itime_n
                E.etimes[n] = curr_itime_n # set the next event time for this node.
            end
        end
    end
end
"""
    epidemic(E,seed,tmax;rseed=-1)

TBW
"""
function epidemic(E,seed,tmax;rseed=-1,debug=false)
    _clear!(E)
    #set random seed
    rseed==-1 ? Random.seed!(abs(rand(Int))) : Random.seed!(rseed)
    _infect_node!(E,0,seed)

    t = 1  
    nnodes = lastindex(E.snodes)
    ninfs = nnodes - sum(E.snodes)
    inf_history = zeros(tmax)
    inf_history[t] = ninfs 
    p = Progress(tmax, dt=1.0)   # minimum update interval: 1 second
    while t<=tmax && ninfs>0 && !isempty(E.etimes) && first(E.etimes)[2]<=tmax
        (n,t) = dequeue_pair!(E.etimes)
        if debug
            println("t=$t")
            println("dequeued ($n,$t)")
            println("total infected nodes: $(nnodes-sum(E.snodes))")
            println("estimated infected nodes: $ninfs")
        end
        if E.snodes[n] && E.itime[n]==t
            _infect_node!(E,t,n)
            ninfs+=1
            inf_history[t]+=1
        else
            _recovery!(E,t,n)
            ninfs-=1
            inf_history[t]-=1
        end
        ProgressMeter.update!(p,t)       
    end
    return inf_history,E.snodes
end

end


using MatrixNetworks
using CairoMakie
using SparseArrays

gname = "study-25-150.smat"
A = MatrixNetworks.readSMAT("data/$gname")

E = EventEpidemic.SISData(A,5e-3,5e-2);
@profview EventEpidemic.epidemic(E,1,100,debug=false);

# fig = Figure()
# ax = Axis(fig[1, 1], 
#     xlabel = "Time Step", 
#     ylabel = "Fraction of Nodes Infected",
#     title="Pairwise-SIS for $(gname[1:end-5])")
# lines!(cumsum(result)./lastindex(A,1))
# ylims!(ax, 0, 1)
# display(fig)



E = EventEpidemic.SISData(A,5e-3,5e-2);
@profview EventEpidemic.epidemic(E,1,100,debug=true);
