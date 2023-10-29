######################################################################################################
# This script is to visualize the arm to reward function μ                                           #
######################################################################################################

using CSV, DataFrames, StatsPlots, Printf;

include("experiment_helpers.jl");

for setting_id=1:5
    μ = generate_settings(setting_id, 40);
    K = length(μ);
    font = Plots.font("Times New Roman", 12)
    fig = plot(collect(1:K), μ, linewidth=1.5, legend=false, xtickfontsize=12, ytickfontsize=12);
    savefig(fig, "setting$(setting_id)_mean.pdf");
end

for setting_id in [6]
    μ = generate_settings(setting_id, 10);
    K = length(μ);
    font = Plots.font("Times New Roman", 12)
    fig = plot(collect(1:K), μ, linewidth=1.5, legend=false, xtickfontsize=12, ytickfontsize=12);
    savefig(fig, "setting$(setting_id)_mean.pdf");
end
