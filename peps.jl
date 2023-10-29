######################################################################################################
# A Pure Exploration problem (pep) is parameterised by
# - a domain, \mathcal M, representing the prior knowledge about the structure
# - a query, as embodied by a correct-answer function istar
# To specify a PE problem here, we need to compute the following things:
# - nanswers: number of possible answers
# - istar: correct answer for feasible μ
######################################################################################################

using IterTools;

"""
Classical.
"""
struct BestArm
    expfam; # common exponential family
end
nanswers(pep::BestArm, μ) = length(μ);
istar(pep::BestArm, μ) = argmax(μ);
getexpfam(pep::BestArm, k) = pep.expfam;
