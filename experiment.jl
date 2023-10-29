using JLD2;
using Distributed;
using Printf;
@everywhere include("runit.jl");
include("experiment_helpers.jl");

setting_id, M = parse(Int, ARGS[1]), parse(Int, ARGS[2]);

# M=K for setting_id∈{1,...,5}, whereas M is as described in Appendix J.6 of the paper
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
        Ts = range(5000,25000,step=5000)
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

for T in Ts
    fname = "setting$(setting_id)_M$(M)_T$(T).dat";
    if !isfile(fname)
        μ = generate_settings(setting_id, M);
        K = length(μ); N = 40000;
        seed = 1234; dist = Bernoulli(); pep = BestArm(dist);
        srs = [
            CR('A'), CR('C'),
            SeqRej(), SeqHalf(),
            UGapEb(),
        ];
        println("μ=$μ,\nK=$K, T=$T, N=$N, seed=$seed");
        @time data = pmap(
            ((sr,i),) -> runit(seed+i, sr, μ, pep, T),
            Iterators.product(srs, 1:N)
        );
        dump_stats(pep, μ, T, srs, data, N);
        @save fname dist μ pep srs data T N seed
    end
end
