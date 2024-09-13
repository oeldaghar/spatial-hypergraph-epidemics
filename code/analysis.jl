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
fnames = get_fnames(get_gname(fnames[1]))
gnames = last.(split.(fnames,'/'))

#plot one - make individual plot for each data set without comparing
fnames = readdir("data/output/sirs/")
@showprogress for fname in fnames 
    data = readdlm(joinpath("data/output/sirs/$fname"))
    parts = split(fname,'%')
    f,ax = make_banded_plot(data,"$(parts[1])\n$(parts[2])")
    plot_file_name = fname[1:end-4]
    plot_file_name="$plot_file_name%banded.png"
    save("data/output/figures/$plot_file_name", f)
end


#plot two - moving alpha 


#plot three - vary normalization term 

# xs = #discrete bins 
# ys = #data 
# boxplot(xs, ys)

function boxplot_alphas_normalization()
    fnames = readdir("data/output/sirs/")
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
        f = boxplot!(ax,Vector{Float64}(xs),Vector{Float64}(ys),width=1/15)
        display(f)
        push!(figures,f)
    end
    return figures 
end


figures = boxplot_alphas_normalization()
display(figures[2])



fnames = readdir("data/output/sirs/")
normalization = "squared"
filtered_fnames = filter(x->occursin("-$normalization",x),fnames)
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
boxplot!(ax,xs,ys,width=1/15)
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
ylims!(ax,100,500)
f


##another one without mixing info across trajectories
f = Figure()
ax = Axis(f[1,1],
        xlabel="Alpha - value of hypergraph interpolate",
        ylabel="trailing infections",
        title="$gname\nNormalization=$normalization")
result = zeros(10,lastindex(alphas))
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
lines!(ax, yavg,color=(:blue,1.0))
band!(ax, alphas, ylow, yhigh, color=(:blue, 0.25))
ylims!(200,450)
f