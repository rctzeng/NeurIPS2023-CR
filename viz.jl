######################################################################################################
# Use this script to plot the experiment results of the compared algorithms for different settings   # 
######################################################################################################
using JLD2, CSV, DataFrames, StatsPlots, Printf;

include("experiment_helpers.jl");
include("runit.jl"); # for types

setting_id, M = parse(Int, ARGS[1]), parse(Int, ARGS[2]);
# M=K for setting_id∈{1,...,6}, whereas M is as described in Appendix J.6 of the paper

df = DataFrame(algorithm=String[], T=Int[], δ=Float32[], K=Int[]);
if setting_id==1
    if M == 10
        Ts = range(4000,8000,step=800)
    elseif M == 20
        Ts = range(6000,16000,step=2000)
    elseif M == 40
        Ts = range(20000,40000,step=5000)
    end
elseif setting_id==2
    if M == 10
        Ts = range(2000,8000,step=1200)
    elseif M == 20
        Ts = range(4000,19000,step=3000)
    elseif M == 40
        Ts = range(5000,25000,step=5000) #range(15000,40000,step=5000)
    end
elseif setting_id==3
    if M == 10
        Ts = range(2000,4000,step=400)
    elseif M == 20
        Ts = range(2000,10000,step=2000)
    elseif M == 40
        Ts = range(15000,25000,step=5000)
    end
elseif setting_id==4
    if M == 10
        Ts = range(400,2000,step=500)
    elseif M == 20
        Ts = range(400,2000,step=500)
    elseif M == 40
        Ts = range(2000,4000,step=400)
    end
elseif setting_id==5
    if M == 10
        Ts = range(500,2500,step=500)
    elseif M == 20
        Ts = range(1500,4000,step=500)
    elseif M == 40
        Ts = range(2000,6000,step=800)
    end
elseif setting_id==6
    if M == 5
        Ts = range(800,2400,step=400)
    elseif M == 6
        Ts = range(1000,2500,step=500)
    elseif M == 10
        Ts = range(1000,5000,step=1000)
    end
end
for Te in Ts
    fname = "setting$(setting_id)_M$(M)_T$(Te).dat";
    if isfile(fname)
        @load fname dist μ pep srs data T N seed
        ⋆ = istar(pep, µ); data = getindex.(data, 1); K = length(μ);
        for r in eachindex(srs)
            es = [maximum([1e-6,err(x[1][1], ⋆)]) for x in data[r,:]];
            δ = mean(es);
            push!(df, (abbrev(srs[r]), T, δ, length(μ)));
        end
    end
end

for T in Ts
    println("[T=$T]");
    for algo in ["CR-A", "CR-C", "SR", "SH", "UG"]
        sdf = filter(row-> (row.T == T) && (row.algorithm == algo), df);
        for delta in sdf.δ
            println("[$algo] δ=$(round(100*delta, digits=3))");
        end
    end
    println("");
end

font = Plots.font("Times New Roman", 12)
y_max, y_min = maximum(df.δ), minimum(df.δ);
y_max, y_min = log10(y_max)+0.1, log10(y_min)-0.1;
fig = plot();
for algo in ["UG", "SH", "SR", "CR-C", "CR-A"]
    sdf = filter(row->row.algorithm == algo, df); sort!(sdf, :T);
    plot!(sdf.T, sdf.δ, xlabel="budget", label=algo, yscale=:log10, xtickfont=font, ytickfont=font, legendfont=font, legend=:outertopright, linewidth=1.5, fillalpha=0.15);
end
ylims!((10^(y_min), 10^(y_max)));
savefig(fig, "setting$(setting_id)_K$(df[1,:].K).pdf");
