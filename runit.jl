######################################################################################################
# Copy from https://bitbucket.org/wmkoolen/tidnabbil/src/master/purex_games_paper/ by Wouter Koolen  #
######################################################################################################
using Random;
using CPUTime;
include("peps.jl");
include("expfam.jl");
include("samplingrules.jl");

# Run the learning algorithm, paramterised by a sampling rule
# The stopping and recommendation rules are common
# Ts must be a list of thresholds *in increasing order*

function runit(seed, sr, μs, pep, T)
    rng = MersenneTwister(seed);
    K = length(μs);
    N = zeros(Int64, K); # counts
    S = zeros(K);        # sum of samples
    baseline = CPUtime_us();
    # pull each arm once
    for k in 1:K
        S[k] += sample(rng, getexpfam(pep, k), μs[k]); N[k] += 1;
    end

    state = start(sr, N, T);
    R = Tuple{Int64, Array{Int64,1}, UInt64}[]; # collect return values
    while true
        t = sum(N);
        hμ = S./N;
        ⋆ = istar(pep, hμ);
        while t >= T
            best_arm = decision(state, pep, ⋆, N, S, T);
            push!(R, (best_arm, copy(N), CPUtime_us()-baseline));
            return R;
        end
        # invoke sampling rule
        k = nextsample(state, pep, ⋆, N, S, T);
        # and actually sample
        reward = sample(rng, getexpfam(pep, k), μs[k]);
        S[k] += reward; N[k] += 1;
    end
end
