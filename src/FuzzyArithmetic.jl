module FuzzyArithmetic

using ProbabilityBoundsAnalysis, IntervalArithmetic, PyPlot

using BivariateCopulas: M, W, Pi, Gaussian

import Base: -, +, *, /, //, <, >, ⊆, ^, intersect, issubset, rand, min, max, log, exp, sin, cos, tan, isequal, ∪, ∩
import ProbabilityBoundsAnalysis: pbox, plot, left, right, mean, var, env, makepbox

abstract type AbstractPoss <: Real end


export FuzzyNumber, Fuzzy, mass, membership, isnested, iscons, mean, cut, descritize, makeCons, env,
    makepbox, makefuzzy, ecdf2fuzzy, isfuzzy, DSS2Fuzzy,

    # arithemtic
    levelwise, levelwiseOpp, supT, mobiusTransform2D, sigmaFuzzy, tauFuzzy,

    # plots
    plot

include("FuzzyNumber.jl")
include("arithmetic.jl")
# include("Tnorms.jl")
# include("PossibilityNumbers.jl")
# include("PossArithmetic.jl")
include("plots.jl")
include("inference.jl")



function __init__()
    using3D()
end

end # module
