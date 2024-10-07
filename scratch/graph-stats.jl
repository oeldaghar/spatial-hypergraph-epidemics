#just checking projected graph stats..

include("../code/hypergraph-tools.jl")

using Arpack
using MatrixNetworks
#uses data loading functiosn from experiments-testing-v1.. move those to separate file and load them
fnames = readdir("data/hypergraphs/")
filter!(x->endswith(x,".txt"),fnames)

gname = get_gname(fnames[1])
fnames = get_fnames(gname)
hedges = read_hedges(fnames[3])


projected_edges = project(hedges)
A = edges2sparse(projected_edges)

nnz(A)/lastindex(A,1)
lastindex(largest_component(A)[1],1)
eigs(A,nev=1)

