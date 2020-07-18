###
#   This file is a part of FuzzyArithmetic.jl package
#
#   This file defines a FuzzyNumber type, and various supporting functions
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
###

using IntervalArithmetic
using ProbabilityBoundsAnalysis
using PyPlot

import Base: -, +, *, /, //, <, >, ⊆, ^, intersect, issubset, rand, min, max

abstract type AbstractFuzzy <: Real end

struct FuzzyNumber <: AbstractFuzzy

    Core :: Interval{<:Real}
    Range :: Interval{<:Real}
    Membership :: Array{Interval{<:Real}, 1}

    function FuzzyNumber(Core = interval(0.5), Range = interval(0, 1), Membership = missing; steps = 200)
        
        if !(Core ⊆ Range); throw(ArgumentError("Core must be a subset of the Range.\nProvided Core = $Core\nProvided Range = $Range")); end

        if ismissing(Membership); Membership = linearInterp(Core, Range, steps); end

        if Range != Membership[1]; throw(ArgumentError("Range must be the first element of the Membership function.\nProvided Range = $Range\nProvided Membership[1] = $(Membership[1]))"));end
        if Core != Membership[end]; throw(ArgumentError("Core must be the last element of the Membership function.\nProvided Core = $Core\nProvided Membership[end] = $(Membership[end]))"));end

        return new(Core, Range, Membership)
    end
end

function (obj::FuzzyNumber)(x, y)
    return mass(obj, x, y)
end


function (obj::FuzzyNumber)(x :: Interval{T}) where T <:Real
    return mass(obj, x)
end

FuzzyNumber(lowerbound :: Real, Core :: Real, upperbound :: Real; steps = 200) = FuzzyNumber(interval(Core), interval(lowerbound, upperbound), steps = steps)
FuzzyNumber(lowerbound :: Real, CoreLeft :: Real, CoreRight :: Real, upperbound :: Real; steps = 200) = FuzzyNumber(interval(CoreLeft, CoreRight), interval(lowerbound, upperbound), steps = steps)
FuzzyNumber(lowerbound :: Real, Core :: Interval{<:Real}, upperbound :: Real; steps = 200) = FuzzyNumber(Core, interval(lowerbound, upperbound), steps = steps)

function FuzzyNumber( Membership :: Array{Interval{T}, 1}) where T <: Real

    if !isnested(Membership); throw(ArgumentError("Invalid Membership function, intervals must be nested from Range (Membership[1]) to Core (Membership[end])")); end
    Core = Membership[end]; Range = Membership[1];
    return FuzzyNumber(Core, Range, Membership, steps = length(Membership))

end


function isnested( a :: Array{Interval{T}, 1}) where T <: Real
    leftSorted = issorted( left.(a) ) 
    rightSorted = issorted( right.(a), rev=true) 
    return leftSorted & rightSorted
end

function mass( x :: FuzzyNumber, lo :: Real , hi :: Real)

    if hi < lo; throw(ArgumentError("Incompatable bounds. Provided: [$lo, $hi]")); end

    if x.Membership[1] ⊂ interval(lo, hi); return interval(1,1); end
    if x.Membership[1].lo > hi; return interval(0,0); end
    if x.Membership[1].hi < lo; return interval(0,0); end

    lefts = left.(x.Membership);
    rights = right.(x.Membership);

    j = range(0, 1, length = length(x.Membership))

    vals = [lefts; reverse(rights)];
    jj = [j;reverse(j)]

    inside = lo .<= vals .<= hi

    Poss = maximum( jj[inside])
    Ness = 1 - maximum(jj[ .~inside])

    return interval(Ness, Poss)

end

function mass( x :: FuzzyNumber, y :: Interval{T}) where T <: Real  
    return mass( x, y.lo, y.hi)
end

function cut( x :: FuzzyNumber, α :: Real )

    if !(α ∈ interval(0, 1) ); throw(ArgumentError("α ∉ [0, 1]\nProvided α = $α")); end

    α = (α *(length( x.Membership ) -1)) + 1; # Scale α by number of elements
    
    i = Int(round(α, RoundDown))         # Round down for conservatism
    return x.Membership[i]
end

#α(x :: FuzzyNumber, α :: Real) = cut(x :: FuzzyNumber, α :: Real)
#alpha(x :: FuzzyNumber, α :: Real) = cut(x :: FuzzyNumber, α :: Real)
#alphaCut(x :: FuzzyNumber, α :: Real) = cut(x :: FuzzyNumber, α :: Real)


function makepbox( x:: FuzzyNumber)
    return pbox()
end


function makefuzzy(x :: pbox)
    return FuzzyNumber()
end

function linearInterp(Core:: Interval, Range :: Interval, steps = 200)
    return interval.(range(Range.lo, Core.lo, length = steps), range(Range.hi, Core.hi, length = steps))
end


function Base.show(io::IO, z::FuzzyNumber)

    Range = z.Range
    Core = z.Core
    if isscalar(z.Core); Core = z.Core.lo; end
    
    print(io, "Fuzzy: \t ~ ( Range=$Range, Core=$Core )");

end

include("plots.jl")
include("FuzzyArithmetic.jl")
