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

#using IntervalArithmetic
using ProbabilityBoundsAnalysis
using PyPlot
using IntervalUnionArithmetic
using3D()
using BivariateCopulas: M, W, Pi, Gaussian

import Base: -, +, *, /, //, <, >, ⊆, ^, intersect, issubset, rand, min, max, log, exp, sin, cos, tan, isequal, ∪
import ProbabilityBoundsAnalysis: pbox, plot, left, right, mean, var, env

abstract type AbstractPoss <: Real end

struct FuzzyNumber <: AbstractPoss

    Core :: Interval{T} where T <: Real
    Range :: Interval{T} where T <: Real
    Membership :: Array{Interval{T}, 1} where T <: Real

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

Fuzzy = FuzzyNumber

FuzzyNumber(lowerbound :: Real, Core :: Real, upperbound :: Real; steps = 200) = FuzzyNumber(interval(Core), interval(lowerbound, upperbound), steps = steps)
FuzzyNumber(lowerbound :: Real, CoreLeft :: Real, CoreRight :: Real, upperbound :: Real; steps = 200) = FuzzyNumber(interval(CoreLeft, CoreRight), interval(lowerbound, upperbound), steps = steps)
FuzzyNumber(lowerbound :: Real, Core :: Interval{<:Real}, upperbound :: Real; steps = 200) = FuzzyNumber(Core, interval(lowerbound, upperbound), steps = steps)

function FuzzyNumber( Membership :: Array{Interval{T}, 1}) where T <: Real

    #if !isnested(Membership); throw(ArgumentError("Invalid Membership function, intervals must be nested from Range (Membership[1]) to Core (Membership[end])")); end
    Core = Membership[end]; Range = Membership[1];
    return FuzzyNumber(Core, Range, Membership, steps = length(Membership))

end

#=
function isnested( a :: Array{Interval{T}, 1}) where T <: Real
    leftSorted = issorted( left.(a)) 
    rightSorted = issorted( right.(a),  rev=true) 
    return leftSorted & rightSorted
end
=#

#isnested(a :: Array{Interval{T}, 1}) where T <: Real = issorted(a, lt=⊂, rev=true)

isnested(a :: Array{Interval{T}, 1}) where T <: Real = all(a[2:end] .⊆ a[1:end-1])
#isnested(a :: Array{Interval{T}, 1}) where T <: Real = true
iscons( a :: Array{Interval{T}, 1}) where T <: Real = isnested(a)

function mass( x :: FuzzyNumber, lo :: Real , hi :: Real)

    if hi < lo; throw(ArgumentError("Incompatable bounds. Provided: [$lo, $hi]")); end

    if x.Membership[1] ⊂ interval(lo, hi); return interval(1,1); end
    if x.Membership[1].lo > hi; return interval(0,0); end
    if x.Membership[1].hi < lo; return interval(0,0); end

    lefts = left.(x.Membership);
    rights = right.(x.Membership);

    j = range(0, 1, length = length(x.Membership)+1)

    j = j[2:end];

    vals = [lefts; reverse(rights)];
    jj = [j;reverse(j)]

    inside = lo .<= vals .<= hi

    Poss = maximum( jj[inside])
    Ness = 1 - maximum(jj[ .~inside])

    return interval(Ness, Poss)

end

function membership(F :: FuzzyNumber, x ::Union{Float64, Int64})
    mems = F.Membership
    here = findlast(x .∈ mems)
    if isnothing(here) return 0;end
    return here/(length(mems)+1)
end

function mass( x :: FuzzyNumber, y :: Interval{T}) where T <: Real  
    return mass( x, y.lo, y.hi)
end

function mean( x :: FuzzyNumber) 
    mems = x.Membership;
    EL = mean(left.(mems))
    EU = mean(right.(mems))
    return interval(EL, EU)
end

function cut( x :: FuzzyNumber, α :: Real )

    if !(α ∈ interval(0, 1) ); throw(ArgumentError("α ∉ [0, 1]\nProvided α = $α")); end

    Mems = deepcopy(x.Membership)
    push!(Mems, Mems[end])

    α = (α *(length( Mems ) -1)) + 1; # Scale α by number of elements
    
    i = Int(round(α, RoundDown))         # Round down for conservatism
    return Mems[i]
end

#α(x :: FuzzyNumber, α :: Real) = cut(x :: FuzzyNumber, α :: Real)
#alpha(x :: FuzzyNumber, α :: Real) = cut(x :: FuzzyNumber, α :: Real)
#alphaCut(x :: FuzzyNumber, α :: Real) = cut(x :: FuzzyNumber, α :: Real)

##
#   Reshapes membership function to match number of steps
##
function descritize(x :: FuzzyNumber, steps = 200)

    is = range(0, 1, length= steps)

    NewMems = cut.(x,is)

    return Fuzzy(NewMems)

end


function makeCons(x :: Array{Interval{T},1}) where T <: Real  

    x = sort(x,lt = ⊂); x = reverse(x)
    z = Interval{Float64}[]
    for i =1:length(x)-1
        lo = x[i+1].lo; hi = x[i+1].hi;
        if x[i].lo < lo; lo = x[i].lo;end
        if x[i].hi > hi; hi = x[i].hi;end
        push!(z, interval(lo, hi))
    end
    push!(z,x[end])
    if !isnested(z) z= makeCons(z);end
    return z
end

function makeConsDom(mem :: Array{Interval{T},1}, masses ) where T <: Real  

    lefts = left.(mem);
    rights = right.(mem);

    subPos = interval(lefts, rights)

    NewMembership = [sum(masses[ mem .⊇ this ]) for this in subPos]  # Find beliefs

    


end

function makeCons1(x :: Array{Interval{T},1}) where T <: Real  

    x = sort(x,lt = ⊂); #x = reverse(x)
    z = Interval{Float64}[]
    for i =1:length(x)-1
        thisZ = x[i]
        if !(x[i] ⊆ x[i+1]); thisZ = hull(x[i],x[i+1]);end
        push!(z, thisZ)
    end
    push!(z,x[end])
    return reverse(z)
end

dia(x, y) = diam(x) < diam(y)

function makeCons2(x :: Array{Interval{T},1}) where T <: Real  

    x = sort(x,rev = true);
    z = Interval{Float64}[]
    for i =1:length(x)-1
        lo = x[i].lo; hi = x[i].hi;
        if x[i].lo > x[i+1].lo; lo = x[i+1].lo;end
        if x[i].hi < x[i+1].hi; hi = x[i+1].hi;end
        push!(z, interval(lo, hi))
    end
    push!(z,x[end])
    return z
end

function env(x :: FuzzyNumber, y :: FuzzyNumber)

    Mems = hull.(x.Membership, y.Membership)
    return FuzzyNumber(Mems)

end


function makepbox( x:: FuzzyNumber)

    lefts = left.(x.Membership)
    rights = right.(x.Membership)

    newRights = reverse(rights)

    return pbox(lefts, newRights)
end


function makefuzzy(x :: pbox)

    # Need to check if possible

    lefts = x.u;
    rights = x.d;

    newRight = reverse(rights);
    
    num = Int(floor(length(lefts)/2));

    Fuzz = FuzzyNumber(interval.(lefts[1:num], newRight[1:num]))

    return descritize(Fuzz,200)
end

function makefuzzy(x :: Union{Real, Interval}, steps = 200) 
    Mems = [interval(x) for i=1:steps]
    return Fuzzy(Mems)
end

function ecdf2fuzzy(x)
    xSort = sort(x)
    l = length(x)
    if iseven(l)
        mid = Integer(l/2)
        lefts = xSort[1:mid]
        rights = reverse(xSort[mid+1:end])
        Mems = interval.(lefts, rights)
    else   
        mid = Integer((l-1)/2)
        lefts = xSort[1:mid]
        rights = reverse(xSort[mid+2:end])
        core = Interval(xSort[mid+1])
        Mems = interval.(lefts,rights)
        push!(Mems, core)
    end
    return Fuzzy(Mems)
end

function linearInterp(Core:: Interval, Range :: Interval, steps = 200)
    return interval.(range(Range.lo, Core.lo, length = steps), range(Range.hi, Core.hi, length = steps))
end


isfuzzy(x) = typeof(x) <: FuzzyNumber

function Base.show(io::IO, z::FuzzyNumber)

    Range = z.Range
    Core = z.Core
    if isscalar(z.Core); Core = z.Core.lo; end
    
    print(io, "Fuzzy: \t ~ ( Range=$Range, Core=$Core )");

end


include("FuzzyArithmetic.jl")
#include("Tnorms.jl")
include("PossibilityNumbers.jl")
include("PossArithmetic.jl")
include("plots.jl")
include("inference.jl")