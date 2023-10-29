using Statistics;

function err(v, ⋆)
    return (v!=⋆) ? 1 : 0;
end

function generate_settings(setting_id, M)
    μ = [];
    if setting_id == 1  # one group of suboptimal arms
        push!(μ, 0.5);
        for k=1:(M-1)
            push!(μ, 0.45);
        end
    elseif setting_id == 2  # two groups of suboptimal arms
        push!(μ, 0.5);
        K1 = floor(Int, (M-1)/2); K2 = M-1-K1;
        for k=1:K1
            push!(μ, 0.45);
        end
        for k=1:K2
            push!(μ, 0.4);
        end
    elseif setting_id == 3  # linear
        push!(μ, 0.75);
        for k=1:(M-1)
            push!(μ, 0.75-0.5*(k/M));
        end
    elseif setting_id == 4  # concave
        push!(μ, sin(π/2*(M-1)/M));
        for k=2:M
            push!(μ, sin(0.45*π*(M-k+1)/M));
        end
    elseif setting_id == 5  # convex
        for k=1:M
            push!(μ, 0.3/(k+1));
        end
    elseif setting_id == 6  # stair
        for m=1:M
            for j=1:m
                push!(μ, 0.75*(1/3)^(m/M));
            end
        end
    elseif setting_id == 7 # 0.99, 0.98
        push!(µ, 0.9);
        push!(µ, 0.8);
        for k=3:M
            push!(µ, 0);
        end
    end
    return μ
end

function dump_stats(pep, μ, T, srs, datas, repeats)
    K = length(μ);
    best_arm = typeof(pep) == BestArm;

    data = getindex.(datas, 1);
    ⋆ = istar(pep, µ);
    rule = repeat("-", 90);
    println("");
    println(rule);
    println("$pep at T = $T");
    println(@sprintf("%23s", "δ"), " ",
            @sprintf("%8s", "σ"), " ",
            @sprintf("%8s", "time"), " ",
            join(map(k -> @sprintf("%15s", k), 1:length(μ))),
    );
    if best_arm
        println(@sprintf("%-33s", "μ"), join(map(x -> @sprintf("     %10.4f", x), μ)));
    end
    println(rule);

    for r in eachindex(srs)
        es = [err(x[1][1], ⋆) for x in data[r,:]];
        size = floor(Int, repeats/5);
        δs = [mean(es[(size*(k-1)+1):(size*k)]) for k=1:5];
        δ, σ = mean(δs), std(δs);
        tim = sum(x->x[3],data[r,:])/repeats;
        println(@sprintf("%-20s", long(srs[r])),
                @sprintf("%0.6f", δ), " ",
                @sprintf("%0.6f", σ), " ",
                @sprintf("%3.1f    ", tim/1e6),
                join(map(k -> @sprintf(" %6.0f(%.4f)",
                    sum(x->x[2][k], data[r,:])/repeats,
                    sum(x->x[2][k], data[r,:])/(repeats*T)), 1:K)),
                " ");
    end
    println(rule);
end
