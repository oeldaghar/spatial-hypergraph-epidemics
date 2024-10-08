#file for analyzing output from experiments
include("data-io.jl")

using CairoMakie
using Statistics
using ProgressMeter

####### INDIVIDUAL PLOTS 
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

function make_individual_figures()
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
end
make_individual_figures()

####### FIGURE 4 
function figure4_hyperedge_sirs_difference()
    #plotting params 
    titlesize = 25
    
    xlabel_size = 25
    xtick_label_size = 20

    ylabel_size = 25
    ytick_label_size = 20

    #load files 
    fnames = readdir("data/output/sirs/")
    filter!(x->endswith(x,".txt"),fnames)
    #filter to the large graphs  
    filter!(x->occursin("-50000-",x),fnames)
    parts = map(x->split(x,"-"),fnames)
    #zoom in on alpha=2.0
    filtered_fnames = filter(x->parse(Float64,split(x,"-")[5])==2.0,fnames)

    f = Figure(resolution=(1400,400))
    ax1 = Axis(f[1,1], 
            title= "g(m)=m",
            xlabel="Time",
            ylabel="Total Infections",
            titlesize = titlesize,
            yticklabelsize = ytick_label_size,
            xticklabelsize = xtick_label_size,
            ylabelsize=ylabel_size,
            xlabelsize=xlabel_size)
    hidespines!(ax1)
    ax2 = Axis(f[1,2], 
                title= "g(m)=sqrt(m)",
                xlabel="Time",
                titlesize = titlesize,
                yticklabelsize = ytick_label_size,
                xticklabelsize = xtick_label_size,
                xlabelsize=xlabel_size,
                yticklabelsvisible=false,
                yticksvisible = false)
    hidespines!(ax2)
    ax3 = Axis(f[1,3], 
            title= "g(m)=m²",
            xlabel="Time",
            titlesize = titlesize,
            yticklabelsize = ytick_label_size,
            xticklabelsize = xtick_label_size,
            xlabelsize=xlabel_size,
            yticklabelsvisible=false,
            yticksvisible = false)
    hidespines!(ax3)
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
    CairoMakie.save(filename,f,px_per_unit=300)
    return f 
end
figure4_hyperedge_sirs_difference()

####### FIGURE 5
function figure5_trailing_infs_plot()
    #plotting params 
    titlesize = 25
    
    xlabel_size = 25
    xtick_label_size = 20

    ylabel_size = 25
    ytick_label_size = 20

    #file names 
    fnames = readdir("data/output/sirs/")
    filter!(x->endswith(x,".txt"),fnames)

    final_fig = Figure(resolution=(1400,400))

    titles = Dict()
    titles["linear"] = "g(m)=m"
    titles["sqrt"] = "g(m)=sqrt(m)"
    titles["squared"] = "g(m)=m²"

    gname = ""

    for (ind,normalization) in enumerate(["linear","sqrt","squared"])
        
        yinfo_bool = ind==1
        ax = Axis(final_fig[1,ind],
            xlabel="α",
            ylabel= yinfo_bool ? "Trailing Infections" : "",
            title=titles[normalization],
            titlesize = titlesize,
            yticklabelsize = ytick_label_size,
            xticklabelsize = xtick_label_size,
            ylabelsize=ylabel_size,
            xlabelsize=xlabel_size,
            yticklabelsvisible = yinfo_bool ? true : false,
            yticksvisible = yinfo_bool ? true : false)
        hidespines!(ax)
        
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
    fname = "data/output/figures/final/trailing-infections-normalization-comparison.pdf"
    CairoMakie.save(fname,final_fig,px_per_unit=300)
    return final_fig
end
figure5_trailing_infs_plot()
