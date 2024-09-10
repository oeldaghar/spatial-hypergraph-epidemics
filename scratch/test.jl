#embed 100 points and add 40 hyperedges with a prescribed radius

using Random
using Makie
using CairoMakie

function plot_circle(ax, center, radius; kwargs...)
    θ = 0:0.01:2π
    x = center[1] .+ radius * cos.(θ)
    y = center[2] .+ radius * sin.(θ)
    lines!(ax, x, y; kwargs...)
end

# Create the scatter plot
function embed_centroids_plot()
    fig = Figure()
    ax = Axis(fig[1, 1])
    scatter!(ax, xy[1,:],xy[2,:])
    scatter!(ax, centroids[1,:],centroids[2,:],color=:black)

    for r in [1.5e-1]
        for x=1:size(centroids)[2]
            c = centroids[:,x]
            plot_circle(ax,c,r,color=:red,alpha=0.5)
        end
    end

    display(fig)
end

d = 2
n = 100000
# k = 20
# r = 5e-2

xy = rand(d,n)
centroids = rand(d,k)

embed_centroids_plot()


fig = Figure()
ax = Axis(fig[1, 1])
scatter!(ax, xy[1,:],xy[2,:])
display(fig)
