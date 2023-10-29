using StatsBase

################################################# Ours #################################################
struct CR
    mode;   # 'A': average, 'C': conservative
end
long(sr::CR) = "CR-$(sr.mode)";
abbrev(sr::CR) = "CR-$(sr.mode)";
mutable struct CRState
    C;      # active set
    mode;   # 'A': average, 'C': conservative
    θ_0;
    CRState(K, mode) = new(collect(1:K), mode, 10^(-5));
end
function start(sr::CR, N, T)
    K = length(N);
    CRState(K, sr.mode);
end
function nextsample(sr::CRState, pep, ⋆, N, S, T)
    K = length(N); 
    if sum(N) < sr.θ_0 * T
        return argmin(N);
    end
    j = length(sr.C); notC = setdiff(collect(1:K), sr.C);
    hμ = [S[k]/N[k] for k=1:K];
    l = sr.C[argmin([hμ[k] for k in sr.C])];
    β = logbar(j)*sum([N[k] for k in sr.C])/(T-sum([N[k] for k in notC]));
    if (j>2) && (isempty(notC) || (!isempty(notC) && N[l] > maximum([N[k] for k in notC])))
        if sr.mode=='A' && (sum([hμ[k] for k in sr.C if k!=l])/(j-1)-hμ[l] > G(β))
            setfield!(sr, :C, filter(v->v!=l, sr.C));
        end
        if sr.mode=='C' && (minimum([hμ[k] for k in sr.C if k!=l])-hμ[l] > G(β))
            setfield!(sr, :C, filter(v->v!=l, sr.C));
        end
    end
    return sr.C[argmin([N[k] for k in sr.C])];
end
function decision(sr::CRState, pep, ⋆, N, S, T)
    return sr.C[argmax([S[k]/N[k] for k in sr.C])];
end
function G(β)
    return 1/sqrt(β)-1;
end
function logbar(j)
    return 0.5 + sum([1/k for k=2:j]);
end

################################################# Baselines #################################################
"""
Sequential Rejection
"""
struct SeqRej end
long(sr::SeqRej) = "SeqRej";
abbrev(sr::SeqRej) = "SR";
mutable struct SeqRejState
    C;  # active set
    p;  # default budget threshold
    SeqRejState(K, p) = new(collect(1:K), p);
end
function start(sr::SeqRej, N, T)
    K = length(N);
    P = 0.5 + sum([1/k for k=2:K]);
    p = [1/(P*k) for k=1:K]; p[1] = p[2];
    SeqRejState(K, p);
end
function nextsample(sr::SeqRejState, pep, ⋆, N, S, T)
    j = length(sr.C);
    CN = [N[k] for k in sr.C];
    n_min = argmin(CN);
    if (CN[n_min] < sr.p[j] * T) || (j==2)
        return sr.C[n_min];
    end
    Cμ = [S[k]/N[k] for k in sr.C];
    l_min = argmin(Cμ);
    setfield!(sr, :C, filter(v->v!=sr.C[l_min], sr.C));
    return nextsample(sr, pep, ⋆, N, S, T);
end
function decision(sr::SeqRejState, pep, ⋆, N, S, T)
    Cμ = [S[k]/N[k] for k in sr.C];
    return sr.C[argmax(Cμ)];
end


"""
Sequential Halving (Karnin et al.)
"""
struct SeqHalf end
long(sr::SeqHalf) = "SeqHalf";
abbrev(sr::SeqHalf) = "SH";
mutable struct SeqHalfState
    C;  # active set
    p;  # budget threshold
    SeqHalfState(K) = new(collect(1:K), [1/(k*log2(K)) for k=1:K]);
end
function start(sr::SeqHalf, N, T)
    SeqHalfState(length(N));
end
function nextsample(sr::SeqHalfState, pep, ⋆, N, S, T)
    j = length(sr.C);
    CN = [N[k] for k in sr.C];
    n_min = argmin(CN);
    if (CN[n_min] < T * sr.p[j]) || (j==2)
        return sr.C[n_min];
    end
    Cμ = [S[k]/N[k] for k in sr.C];
    Cidx = sortperm(Cμ);
    setfield!(sr, :C, [sr.C[Cidx[j-i+1]] for i=1:ceil(Int, j/2)]);
    return nextsample(sr, pep, ⋆, N, S, T)
end
function decision(sr::SeqHalfState, pep, ⋆, N, S, T)
    Cμ = [S[k]/N[k] for k in sr.C];
    return sr.C[argmax(Cμ)];
end


"""
UGapEb (Gabillon et al.)
"""
struct UGapEb end
long(sr::UGapEb) = "UGapEb";
abbrev(sr::UGapEb) = "UG";
struct UGapEbState end
function start(sr::UGapEb, N, T)
    return sr;
end
function confidence_bound(a, b, N, hμ, K)
    β = [b*sqrt(a/N[k]) for k=1:K];
    U = [hμ[k]+β[k] for k=1:K];
    L = [hμ[k]-β[k] for k=1:K];
    B = [maximum([U[j] for j=1:K if j!=k])-L[k] for k=1:K];
    return β, U, L, B
end
function nextsample(sr::UGapEb, pep, ⋆, N, S, T) # sampling rule
    K = length(N);
    # compute the complexity term H by empirical mean (Audibert and Bubeck 2010)
    hμ = S./N; H = sum([4/(hμ[⋆]-hμ[k])^2 for k=1:K if k!=⋆]); a = (T-K)/(4*H);
    β, U, L, B = confidence_bound(a, 1, N, hμ, K);
    J = argmin(B); notJ = [k for k=1:K if k!=J];
    u = notJ[argmax([U[k] for k in notJ])];
    return (β[J]>β[u]) ? J : u;
end
function decision(sr::UGapEb, pep, ⋆, N, S, T) # decision rule
    K = length(N); hμ = S./N;
    H = sum([4/(hμ[⋆]-hμ[k])^2 for k=1:K if k!=⋆]); a = (T-K)/(4*H);
    β, U, L, B = confidence_bound(a, 1, N, S./N, length(N));
    return argmin(B);
end


"""
Uniform sampling
"""
struct RoundRobin end
long(sr::RoundRobin) = "Uniform";
abbrev(sr::RoundRobin) = "RR";
function start(sr::RoundRobin, N, T)
    return sr;
end
function nextsample(sr::RoundRobin, pep, ⋆, N, S, T)
    return 1+(sum(N) % length(N));
end
function decision(sr::RoundRobin, pep, ⋆, N, S, T)
    return ⋆;
end
