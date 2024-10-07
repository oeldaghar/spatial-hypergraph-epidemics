#file for analyzing output from experiments
include("data-io.jl")

using CairoMakie
using Statistics
using ProgressMeter

function make_banded_plot(data,title="")
    # Make plotting data 
    f = Figure()
    ax = Axis(f[1,1], 
                title=title,
                xlabel="Time",
                ylabel="Total Infections")

    hyper_avgs = mapslices(x->mean(x),data,dims=1)
    hyper_ylow = mapslices(x->quantile(x,0),data,dims=1)
    hyper_yhigh = mapslices(x->quantile(x,1),data,dims=1)
    hyper_avgs,hyper_ylow,hyper_yhigh = vec(hyper_avgs),vec(hyper_ylow),vec(hyper_yhigh)

    #render 
    lines!(ax, hyper_avgs,color=(:blue,1.0))
    band!(ax, 1:lastindex(hyper_avgs), hyper_ylow, hyper_yhigh, color=(:blue, 0.25))
    return f,ax 
end

fnames = filter!(x->endswith(x,".txt"),readdir("data/hypergraphs/"))
fnames = get_fnames(get_gname(fnames[end]))
filter!(x->endswith(x,".txt"),fnames)
gnames = last.(split.(fnames,'/'))

#plot one - make individual plot for each data set without comparing
fnames = readdir("data/output/sirs/")
filter!(x->endswith(x,".txt"),fnames)
@showprogress for fname in fnames 
    data = readdlm(joinpath("data/output/sirs/$fname"))
    parts = split(fname,'%')
    f,ax = make_banded_plot(data,"$(parts[1])\n$(parts[2])")
    plot_file_name = fname[1:end-4]
    plot_file_name="$plot_file_name%banded.png"
    save("data/output/figures/$plot_file_name", f)
end

#zoom in on alpha = 2 to highlight hypergraph effects 
fnames = readdir("data/output/sirs/")
filter!(x->endswith(x,".txt"),fnames)
filter!(x->occursin("-50000-",x),fnames)
parts = map(x->split(x,"-"),fnames)
filtered_fnames = filter(x->parse(Float64,split(x,"-")[5])==2.0,fnames)



f = Figure(resolution=(1400,400))
ax1 = Axis(f[1,1], 
        xlabel="Time",
        ylabel="Total Infections")
ax2 = Axis(f[1,2], 
            xlabel="Time",
            ylabel="Total Infections")
ax3 = Axis(f[1,3], 
        xlabel="Time",
        ylabel="Total Infections")
linkaxes!(ax1,ax2,ax3)

fname = filtered_fnames[1]
data = readdlm(joinpath("data/output/sirs/$fname"))
parts = split(fname,'%')

hyper_avgs = mapslices(x->mean(x),data,dims=1)
hyper_ylow = mapslices(x->quantile(x,0),data,dims=1)
hyper_yhigh = mapslices(x->quantile(x,1),data,dims=1)
hyper_avgs,hyper_ylow,hyper_yhigh = vec(hyper_avgs),vec(hyper_ylow),vec(hyper_yhigh)
#render 
lines!(ax1, hyper_avgs,color=(:blue,1.0))
band!(ax1, 1:lastindex(hyper_avgs), hyper_ylow, hyper_yhigh, color=(:blue, 0.25))    

#second one 
fname = filtered_fnames[2]
data = readdlm(joinpath("data/output/sirs/$fname"))
parts = split(fname,'%')

hyper_avgs = mapslices(x->mean(x),data,dims=1)
hyper_ylow = mapslices(x->quantile(x,0),data,dims=1)
hyper_yhigh = mapslices(x->quantile(x,1),data,dims=1)
hyper_avgs,hyper_ylow,hyper_yhigh = vec(hyper_avgs),vec(hyper_ylow),vec(hyper_yhigh)
#render 
lines!(ax2, hyper_avgs,color=(:blue,1.0))
band!(ax2, 1:lastindex(hyper_avgs), hyper_ylow, hyper_yhigh, color=(:blue, 0.25))    

#third one 
fname = filtered_fnames[3]
data = readdlm(joinpath("data/output/sirs/$fname"))
parts = split(fname,'%')

hyper_avgs = mapslices(x->mean(x),data,dims=1)
hyper_ylow = mapslices(x->quantile(x,0),data,dims=1)
hyper_yhigh = mapslices(x->quantile(x,1),data,dims=1)
hyper_avgs,hyper_ylow,hyper_yhigh = vec(hyper_avgs),vec(hyper_ylow),vec(hyper_yhigh)
#render 
lines!(ax3, hyper_avgs,color=(:blue,1.0))
band!(ax3, 1:lastindex(hyper_avgs), hyper_ylow, hyper_yhigh, color=(:blue, 0.25))    
CairoMakie.xlims!(1,620)
f
filename = "data/output/figures/final/time-average-diffuison-50000-alpha-2.pdf"
CairoMakie.save(filename,f)
# Make plotting data 

###########

for (i,fname) in enumerate(filtered_fnames)
    data = readdlm(joinpath("data/output/sirs/$fname"))
    parts = split(fname,'%')
    f,ax = make_banded_plot(data,"$(parts[1])\n$(parts[2])")
end





#pick specific graph as example 
#linear
#sqrt 
#squared 

fnames
#plot two - moving alpha 


#plot three - vary normalization term 

# xs = #discrete bins 
# ys = #data 
# boxplot(xs, ys)

function boxplot_alphas_normalization()
    fnames = readdir("data/output/sirs/")
    filter!(x->occursin("-50000-",x),fnames)
    figures = []
    for normalization in ["linear","squared","sqrt"]
        filtered_fnames = filter(x->occursin("-$normalization",x),fnames)
        alphas = map(x->x[5],split.(first.(split.(filtered_fnames,'%')),'-'))
        p = sortperm(alphas,by=x->parse(Float64,x))
        filtered_fnames = filtered_fnames[p]
        alphas = alphas[p]
        
        tmp = split.(first.(split.(filtered_fnames,'%')),"-")
        gname = unique(map(x->join(x[1:4],"-"),tmp))

        f = Figure()
        ax = Axis(f[1,1],
                xlabel="Alpha - value of hypergraph interpolate",
                ylabel="trailing average infections",
                title="$gname\nNormalization=$normalization")
        xs = []
        ys = []
        for (ind,fname) in enumerate(filtered_fnames)
            data = readdlm(joinpath("data/output/sirs/$fname"))
            tmp = data[end-1000:end]
            ys = vcat(ys,tmp)
            xs = vcat(xs,[parse(Float64,alphas[ind]) for i=1:lastindex(tmp)])
        end
        f = CairoMakie.boxplot!(ax,Vector{Float64}(xs),Vector{Float64}(ys),width=1/15)
        display(f)
        push!(figures,f)
    end
    return figures 
end


figures = boxplot_alphas_normalization()
display(figures[2])



fnames = readdir("data/output/sirs/")
filter!(x->endswith(x,".txt"),fnames)
normalization = "squared"
filtered_fnames = filter(x->occursin("-$normalization",x),fnames)
filtered_fnames = filter(x->occursin("-50000-",x),filtered_fnames)
alphas = map(x->x[5],split.(first.(split.(filtered_fnames,'%')),'-'))
p = sortperm(alphas,by=x->parse(Float64,x))
filtered_fnames = filtered_fnames[p]
alphas = alphas[p]

tmp = split.(first.(split.(filtered_fnames,'%')),"-")
gname = unique(map(x->join(x[1:4],"-"),tmp))[1]

f = Figure()
ax = Axis(f[1,1],
        xlabel="Alpha - value of hypergraph interpolate",
        ylabel="trailing infections",
        title="$gname\nNormalization=$normalization")
xs = Vector{Float64}()
ys = Vector{Float64}()
for (ind,fname) in enumerate(filtered_fnames)
    data = readdlm(joinpath("data/output/sirs/$fname"))
    tmp = vec(data[1:end,end-1000:end])
    ys = vcat(ys,tmp)
    xs = vcat(xs,[parse(Float64,alphas[ind]) for i=1:lastindex(tmp)])
end
CairoMakie.boxplot!(ax,xs,ys,width=1/15)
# ylims!(ax,100,800)
f


f = Figure()
ax = Axis(f[1,1],
        xlabel="Alpha - value of hypergraph interpolate",
        ylabel="trailing infections",
        title="$gname\nNormalization=$normalization")
xs = Vector{Float64}()
ys = Vector{Float64}()
for (ind,fname) in enumerate(filtered_fnames)
    data = readdlm(joinpath("data/output/sirs/$fname"))
    tmp = data[1:end,end-1000:end]
    ys = push!(ys,mean(tmp))
    xs = push!(xs,parse(Float64,alphas[ind]))
end
lines!(ax,xs,ys)
# ylims!(ax,100,500)
f


fnames = readdir("data/output/sirs/")
filter!(x->endswith(x,".txt"),fnames)


##another one without mixing info across trajectories
function make_trailing_infs_plot(test::Bool=false)
    fnames = readdir("data/output/sirs/")
    filter!(x->endswith(x,".txt"),fnames)

    final_fig = Figure(resolution=(1400,400))
    figures = []
    for (ind,normalization) in enumerate(["linear","sqrt","squared"])
        ax = Axis(final_fig[1,ind],
            xlabel="α",
            ylabel="Average Trailing Infections",
            title="",
            titlesize=24,
            xlabelsize=20,
            ylabelsize=20)
        
        #handle filenames and parameters 
        filtered_fnames = filter(x->occursin("-$normalization",x),fnames)
        filtered_fnames = filter(x->occursin("-50000-",x),filtered_fnames)
        alphas = map(x->x[5],split.(first.(split.(filtered_fnames,'%')),'-'))
        p = sortperm(alphas,by=x->parse(Float64,x))
        filtered_fnames = filtered_fnames[p]
        alphas = alphas[p]

        tmp = split.(first.(split.(filtered_fnames,'%')),"-")
        gname = unique(map(x->join(x[1:4],"-"),tmp))[1]

        result = zeros(25,lastindex(alphas))
        for (ind,fname) in enumerate(filtered_fnames)
            data = readdlm(joinpath("data/output/sirs/$fname"))
            tmp = data[1:end,end-1000:end]
            avgs = vec(mapslices(x->mean(x),tmp,dims=2))
            result[:,ind] = avgs
        end

        yavg = vec(mapslices(x->mean(x),result,dims=1))
        ylow = vec(mapslices(x->quantile(x,0),result,dims=1))
        yhigh = vec(mapslices(x->quantile(x,1),result,dims=1))

        #render 
        lines!(ax, parse.(Float64,alphas), yavg,color=(:blue,1.0),linestyle=:dot)
        CairoMakie.scatter!(ax, parse.(Float64,alphas), yavg,color=(:blue,1.0))
        band!(ax, parse.(Float64,alphas), ylow, yhigh, color=(:blue, 0.25))
        CairoMakie.ylims!(2500,4500)
    end
    return final_fig
end


figs = make_trailing_infs_plot(false)

figs
CairoMakie.save("data/output/figures/final/$gname-normalization-comparison.pdf",figs)
CairoMakie.save("data/output/figures/final/$gname-normalization-comparison.png",figs,px_per_unit = 300)


figs[1]
figs[2]
figs[3]



## full banded simulation 


figs = make_infs_plot()
figs

fnames = readdir("data/output/sirs/")
filter!(x->endswith(x,".txt"),fnames)

figures = []
for normalization in ["linear","squared","sqrt"]
    filtered_fnames = filter(x->occursin("-$normalization",x),fnames)
    filter!(x->occursin(x,"-50000-"),filtered_fnames)
    
    alphas = map(x->x[5],split.(first.(split.(filtered_fnames,'%')),'-'))
    p = sortperm(alphas,by=x->parse(Float64,x))
    filtered_fnames = filtered_fnames[p]
    alphas = alphas[p]
    
    tmp = split.(first.(split.(filtered_fnames,'%')),"-")
    gname = unique(map(x->join(x[1:4],"-"),tmp))

    f = Figure()
    ax = Axis(f[1,1],
            xlabel="Alpha - value of hypergraph interpolate",
            ylabel="trailing average infections",
            title="$gname\nNormalization=$normalization")
    xs = Vector{Float64}()
    ys = Vector{Float64}()
    for (ind,fname) in enumerate(filtered_fnames)
        data = readdlm(joinpath("data/output/sirs/$fname"))
        tmp = data[end-1000:end]
        ys = vcat(ys,tmp)
        xs = vcat(xs,[parse(Float64,alphas[ind]) for i=1:lastindex(tmp)])
    end
    f = CairoMakie.boxplot!(ax,xs,ys,width=1/25)
    push!(figures,f)
end






using Plots

# Sample data
x = 1:10
y_mean = rand(10) * 10  # Mean values
y_lower = y_mean .- 2   # Lower bound (mean - 2)
y_upper = y_mean .+ 2   # Upper bound (mean + 2)

# Create a band plot
Plots.band(x, y_lower, y_upper, label="Confidence Interval", alpha=0.3, color=:blue)

# Optional: Add a line in the middle for the mean
plot!(x, y_mean, label="Mean", color=:black, linewidth=2)

# Add titles and labels
title!("Band Plot Example")
xlabel!("X-axis")
ylabel!("Y-axis")

# Display the plot
display(plot!)
