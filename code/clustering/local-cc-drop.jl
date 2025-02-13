using Random, Plots, LaTeXStrings, Colors
using Measures
using LazySets
using LinearAlgebra

# circle 1
ang1 = range(-pi/4,pi/4,10)
r1 = 1.0
p1 = [0.0;0]
pts1 = [p1[1] .+ r1.*cos.(ang1'); p1[2] .+ r1.*sin.(ang1')]

Plots.scatter([p1[1]],[p1[2]])
Plots.scatter!(pts1[1,:],pts1[2,:])

# circle 2
local b = [Ball2(X[:, v], 0.0005) for v in edge]
local c = ConvexHullArray(b)
Plots.plot!(plt, c, 1e-3, alpha=0.05, c=:blue)