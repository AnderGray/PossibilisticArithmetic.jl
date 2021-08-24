###
#   This file is a part of PossibilisticArithmetic.jl package for performing
#   rigorous operations with possibility distributions/ fuzzy numbers
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
###


module PossibilisticArithmetic

using ProbabilityBoundsAnalysis, IntervalArithmetic, PyPlot

using BivariateCopulas: M, W, Pi, Gaussian

import Base: -, +, *, /, //, <, >, ⊆, ^, intersect, issubset, rand, min, max, log, exp, sin, cos, tan, isequal, ∪, ∩
import ProbabilityBoundsAnalysis: pbox, plot, left, right, mean, var, env, makepbox

abstract type AbstractPoss <: Real end

export FuzzyNumber, Fuzzy, mass, membership, isnested, iscons, mean, cut, descritize, makeCons, env,
    makepbox, makefuzzy, ecdf2fuzzy, isfuzzy, DSS2Fuzzy, check_inside,

    # arithemtic
    levelwise, levelwiseOpp, supT, mobiusTransform2D, sigmaFuzzy, tauFuzzy,

    # plots
    plot,

    left, right

    # Copulas
    M, W, Pi, Gaussian

include("FuzzyNumber.jl")
include("arithmetic.jl")
# include("Tnorms.jl")
# include("PossibilityNumbers.jl")
# include("PossArithmetic.jl")
include("plots.jl")
#include("inference.jl")

left(x :: Fuzzy) = left.(x.Membership)
right(x :: Fuzzy) = right.(x.Membership)

function __init__()
    using3D()
end

end # module
