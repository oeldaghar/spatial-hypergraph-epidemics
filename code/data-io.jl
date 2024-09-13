#TODO file/data naming/loading to a data IO file.
using DelimitedFiles
function get_gname(fname::String)
    parts = split(fname,'-')
    gname = join(parts[1:end-2],'-')
    return gname
end
function get_fnames(gname::String)
    fnames = filter(x->endswith(x,".txt"),readdir("data/hypergraphs/"))

    #filter down to gnames 
    gnames = filter(x->startswith(x,gname),fnames)
    #sort by alpha values 
    sort1 = sort(gnames,by=x->parse(Float64,split(x,'-')[end-1]))
    sort2 = sort(gnames,by=x->parse(Float64,split(x,'-')[5]))

    #checking if file names change 
    if !(sort1==sort2)
        println("WARNING: get_fnames broken, sorting by alpha incorrect, possible change in file name schema or duplicate named graph")
    end
    return ["data/hypergraphs/$x" for x in sort1]
end
function _get_alpha(gname)
    parts = split(gname,'-')
    alpha_v1 = parse(Float64,parts[5])
    alpha_v2 = parse(Float64,parts[end-1])
    if alpha_v1!=alpha_v2
        println("WARNING: _get_alpha broken, possible change in file name schema.")
    end
    return alpha_v1
end
#example usage 
# fnames = filter!(x->endswith(x,".txt"),readdir("data/hypergraphs/"))
# fnames = get_fnames(get_gname(fnames[1]))


# #functions for loading output data from sirs experiments. 
# #fnames look something like "{gname}-{parameters}%{sirs}-{parameters}.txt with {%} serving as the separator between the two
# function get_data_fnames(gname,hyper_beta_func)
#     #get all filenames associated with a graph for a given hyper_beta_funciton (linear, sqrt, squared)
#     fnames = readdir("data/output/sirs")
#     filter!(x->startswith(x,gname),fnames)
#     parts = split.(fnames,'%')
#     gnames = first(parts)
#     epidemics = map(x->x[2],parts)


# end

# tmp = readdir("data/output/sirs")
# g1 = split.(tmp,'%')[1][1]
# filter!(x->startswith(x,g1),tmp)

# println(tmp[1])