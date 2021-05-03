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

"""
    struct FuzzyNumber <: AbstractPoss

Basic type in FuzzyArithmetic. Defines a set of nested intervals (Membership), starting from Range (Membership[1]) to the core (Membership[end])

Used to bound a set of probability distribution functions in the Imprecise Probability sense

# Constructors
* `FuzzyNumber(Core :: Interval, Range :: Interval; steps :: Integer)                   => Give Core and Range, Membership is a linear interpolation`
* `FuzzyNumber(Membership :: Array{Interval{T}, 1} where T <: Real                      => Give Membership`
* `FuzzyNumber(lowerbound :: Real, Core :: Real, upperbound :: Real; steps :: Integer)      => Give Range as Reals`
* `FuzzyNumber(lowerbound :: Real, Core :: Real, upperbound :: Real; steps :: Integer)      => Give Core and Range as reals`

Alias: Fuzzy

See also: [`mass`](@ref), [`membership`](@ref), [`Fuzzy`](@ref), [`plot`](@ref), [`cut`](@ref), [`mean`](@ref)

"""
struct FuzzyNumber <: AbstractPoss

    Membership :: Array{Interval{T}, 1} where T <: Real

    function FuzzyNumber(Core = interval(0.5), Range = interval(0, 1), Membership = missing; steps = 200)

        if ismissing(Membership); Membership = linearInterp(Core, Range, steps+1)[1:steps]; end

        return new(Membership)
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


function mass( x :: FuzzyNumber, y :: Interval{T}) where T <: Real
    Mems = x.Membership;

    masses = 1/interval(length(Mems))

    subs = sum(y .⊂ Mems)
    intersects = sum( y .∩ Mems .!= ∅)

    prob = interval(masses * intersects, masses * subs)

    prob = max(prob, 0)
    prob = min(prob, 1)

    return prob

end

mass( x :: FuzzyNumber, lo :: Real , hi :: Real) = mass(x, interval(lo, hi))

function membership(F :: FuzzyNumber, x ::Union{Float64, Int64})
    mems = F.Membership
    here = findlast(x .∈ mems)
    if isnothing(here) return 0;end
    return here/(length(mems))
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


function makeConsOld(x :: Array{Interval{T},1}) where T <: Real

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


##
#   Imprecise probability to possibility transform
##
function DSS2Fuzzy(FE :: Array{Interval{T},1}, masses; steps = 200) where T <: Real

    if length(FE) < steps; steps = length(FE); end

    lefts = FuzzyArithmetic.left.(FE);
    rights = FuzzyArithmetic.right.(FE);

    lefts = sort(lefts); rights = sort(rights, rev=true)

    subPos = interval.(lefts, rights)

    subPos = unique(subPos)

    NewMembership = [sum(masses[ FE .⊇ this ]) for this in subPos]  # Find beliefs

    memberships = range(0,1,length = steps+1)[1:end-1]

    these = [findall( α .<= NewMembership) for α in memberships]

    return Fuzzy([hull(subPos[this]) for this in these])
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

    Range = z.Membership[1]
    Core = z.Membership[end]
    if isscalar(Core); Core = Core.lo; end

    print(io, "Fuzzy: \t ~ ( Range=$Range, Core=$Core )");

end
