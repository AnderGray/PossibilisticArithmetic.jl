###
#   This file is a part of PossibilisticArithmetic.jl package
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

Basic type in PossibilisticArithmetic. Defines a set of nested intervals (Membership), starting from Range (Membership[1]) to the core (Membership[end])

Used to bound a set of probability distribution functions in the Imprecise Probability sense

# Constructors
* `FuzzyNumber(Core :: Interval, Range :: Interval; steps :: Integer)                   => Give Core and Range, Membership is a linear interpolation`
* `FuzzyNumber(Membership :: Vector{Interval{T}} where T <: Real                        => Give Membership`
* `FuzzyNumber(lowerbound :: Real, Core :: Real, upperbound :: Real; steps :: Integer)      => Give Range as Reals`
* `FuzzyNumber(lowerbound :: Real, Core :: Real, upperbound :: Real; steps :: Integer)      => Give Core and Range as reals`

Alias: Fuzzy

See also: [`mass`](@ref), [`membership`](@ref), [`Fuzzy`](@ref), [`plot`](@ref), [`cut`](@ref), [`mean`](@ref)

"""
struct FuzzyNumber <: AbstractPoss

    Membership::Vector{Interval{T}} where T <: Real

    function FuzzyNumber(Core = interval(0.5), Range = interval(0, 1), Membership = missing; steps = 200)

        if ismissing(Membership); Membership = linearInterp(Core, Range, steps + 1)[1:steps]; end

        return new(Membership)
    end
end

function (obj::FuzzyNumber)(x, y)
    return mass(obj, x, y)
end

function (obj::FuzzyNumber)(x::Interval{T}) where T <: Real
    return mass(obj, x)
end


Fuzzy = FuzzyNumber

FuzzyNumber(lowerbound::Real, Core::Real, upperbound::Real; steps = 200) = FuzzyNumber(interval(Core), interval(lowerbound, upperbound), steps = steps)
FuzzyNumber(lowerbound::Real, CoreLeft::Real, CoreRight::Real, upperbound::Real; steps = 200) = FuzzyNumber(interval(CoreLeft, CoreRight), interval(lowerbound, upperbound), steps = steps)
FuzzyNumber(lowerbound::Real, Core::Interval{<:Real}, upperbound::Real; steps = 200) = FuzzyNumber(Core, interval(lowerbound, upperbound), steps = steps)

function FuzzyNumber(Membership::Vector{Interval{T}}) where T <: Real

    # if !isnested(Membership); throw(ArgumentError("Invalid Membership function, intervals must be nested from Range (Membership[1]) to Core (Membership[end])")); end
    Core = Membership[end]; Range = Membership[1];
    return FuzzyNumber(Core, Range, Membership, steps = length(Membership))

end

#=
function isnested( a :: Array{Interval{T}, 1}) where T <: Real
    leftSorted = issorted( left.(a))
    rightSorted = issorted( right.(a),  rev=true)
    return leftSorted & rightSorted
end =#

# isnested(a :: Vector{Interval{T}}) where T <: Real = issorted(a, lt=⊂, rev=true)

isnested(a::Vector{Interval{T}}) where T <: Real = all(a[2:end] .⊆ a[1:end - 1])
iscons( a::Vector{Interval{T}}) where T <: Real = isnested(a)


function mass(x::FuzzyNumber, y::Interval{T}) where T <: Real
    Mems = x.Membership;

    masses = 1 / interval(length(Mems))

    subs = sum(Mems .⊆ y)
    intersects = sum(y .∩ Mems .!= ∅)

    prob = masses * interval( subs, intersects)

    prob = max(prob, 0)
    prob = min(prob, 1)

    return prob

end

mass( x::FuzzyNumber, lo::Real , hi::Real) = mass(x, interval(lo, hi))

function membership(F::FuzzyNumber, x::Union{Float64,Int64})
    mems = F.Membership
    here = findlast(x .∈ mems)
    if isnothing(here) return 0;end
    return here / (length(mems))
end

function mean(x::FuzzyNumber)
    mems = x.Membership;
    EL = mean(left.(mems))
    EU = mean(right.(mems))
    return interval(EL, EU)
end

function cut(x::FuzzyNumber, α::Real)

    if !(α ∈ interval(0, 1) ); throw(ArgumentError("α ∉ [0, 1]\nProvided α = $α")); end

    Mems = deepcopy(x.Membership)
    push!(Mems, Mems[end])

    α = (α * (length(Mems) - 1)) + 1; # Scale α by number of elements

    i = Int(round(α, RoundDown))         # Round down for conservatism
    return Mems[i]
end

# α(x :: FuzzyNumber, α :: Real) = cut(x :: FuzzyNumber, α :: Real)
# alpha(x :: FuzzyNumber, α :: Real) = cut(x :: FuzzyNumber, α :: Real)
# alphaCut(x :: FuzzyNumber, α :: Real) = cut(x :: FuzzyNumber, α :: Real)

##
#   Reshapes membership function to match number of steps
##
function descritize(x::FuzzyNumber, steps = 200)

    is = range(0, 1, length = steps)

    NewMems = cut.(x, is)

    return Fuzzy(NewMems)

end

function isThisEmpty(x :: Interval)
    if x == ∅ return true; end
    return false
end


function isThisEmpty(x :: IntervalUnion)
    if all(x.v .== ∅) return true; end
    return false
end


##
#   Imprecise probability to possibility transform. Converts a general random set to a fuzzy number.
##

function DSS2Fuzzy(FE::Vector{Interval{T}}, masses::Vector{Float64} = ones(Integer(length(FE))) ./ length(FE); steps = length(FE)) where T <: Real

    if length(FE) < steps; steps = length(FE); end

    lefts = PossibilisticArithmetic.left.(FE);
    rights = PossibilisticArithmetic.right.(FE);

    lefts = sort(lefts); rights = sort(rights, rev = true)   ## Problem ... Would need to mult by 2. But will be over-conservative.

    needReverse = rights .< lefts

    tmp = rights[needReverse]

    rights[needReverse] .= lefts[needReverse]
    lefts[needReverse] .= tmp

    subPos = unique(subPos)
    #subPosSmaller = interval.(nextfloat.(lefts), prevfloat.(rights))    # Required due to rounding in interval arithemtic

    subPos = unique(subPos)
    #subPosSmaller = unique(subPosSmaller)

    NewMembership = [sum(masses[ .!(FE .⊂ this) ]) for this in subPos]  # Find beliefs
    #NewMembership = [sum(masses[ isThisEmpty.(this .∩ FE)]) for this in subPos]  # Find beliefs
    #NewMembership = reverse(NewMembeplotrship)
    zPos = range(0, 1, length=steps+1)

    zMems = Interval{Float64}[]

    push!(zMems, subPos[1])

    for i = 2:length(zPos)-1                                # Find alpha cuts for return Fuzzy
        theseOnes =  zPos[i] .<= NewMembership .< zPos[i+1]
        if !any(theseOnes);
            thisInt = zMems[end]
        else
            thisInt = hull(subPos[theseOnes])
        end
        push!(zMems, thisInt);
    end

    return Fuzzy(zMems)
end

function DSS2FuzzySlow(FE::Vector{Interval{T}}, masses::Vector{Float64} = ones(Integer(length(FE))) ./ length(FE); steps = length(FE)) where T <: Real

    if length(FE) < steps; steps = length(FE); end

    lefts = PossibilisticArithmetic.left.(FE);
    rights = PossibilisticArithmetic.right.(FE);

    lefts = sort(lefts); rights = sort(rights, rev = true)   ## Problem ... Would need to mult by 2. But will be over-conservative.

    subPos = interval.(lefts, rights)

    subPos = unique(subPos)

    subPosC = IntervalUnionArithmetic.complement.(subPos)

    NewMembership = [sum(masses[ .!isThisEmpty.(this .∩ FE) ]) for this in subPosC]  # Find beliefs

    zPos = range(0,1, length=steps+1)

    zMems = Interval{Float64}[]

    push!(zMems, subPos[1])

    for i = 2:length(zPos)-1                               # Find alpha cuts for return Fuzzy
        theseOnes =  zPos[i] .<= NewMembership .< zPos[i+1]
        if !any(theseOnes);
            thisInt = zMems[end]
        else
            thisInt = hull(subPos[theseOnes])
        end
        push!(zMems, thisInt);
    end

    return Fuzzy(zMems)
end


function DSS2Fuzzy2(FE::Vector{Interval{T}}, masses::Vector{Float64} = ones(Integer(length(FE))) ./ length(FE); steps = length(FE)) where T <: Real

    if length(FE) < steps; steps = length(FE); end

    lefts = PossibilisticArithmetic.left.(FE);
    rights = PossibilisticArithmetic.right.(FE);

    lefts = sort(lefts); rights = sort(rights, rev = true)   ## Problem ... Would need to mult by 2. But will be over-conservative.

    subPos = interval.(lefts, rights)

    subPos = unique(subPos)

    NewMembership = [sum(masses[ this .⊆ FE ]) for this in subPos]  # Find beliefs
    zPos = range(0,1, length=steps+1)

    zMems = Interval{Float64}[]

    push!(zMems, subPos[1])

    for i = 2:length(zPos)-1                                # Find alpha cuts for return Fuzzy
        theseOnes =  zPos[i] .<= NewMembership .< zPos[i+1]
        if !any(theseOnes);
            thisInt = zMems[end]
        else
            thisInt = hull(subPos[theseOnes])
        end
        push!(zMems, thisInt);
    end

    return Fuzzy(zMems)
end



function DSS2FuzzyDom(FE::Vector{Interval{T}}, masses::Vector{Float64} = ones(Integer(length(FE))) ./ length(FE); steps = length(FE)) where T <: Real

    if length(FE) < steps; steps = length(FE); end


    widths = diam.(FE);

    idns = sortperm(widths);

    #core = FE[idns[1]]

    lefts = PossibilisticArithmetic.left.(FE);
    rights = PossibilisticArithmetic.right.(FE);

    lCon = lefts[idns];
    rCon = rights[idns];

    Ns = length(FE)
    lConsNew = zeros(Ns)
    rConsNew = zeros(Ns)

    lConsNew[1] = lCon[1]
    rConsNew[1] = rCon[1]

    for i = 2:Ns
        lConsNew[i] = min(lCon[i],lConsNew[i-1])
        rConsNew[i] = max(rCon[i],rConsNew[i-1])

    end



    subPos = interval.(lConsNew, rConsNew)

    subPos = unique(subPos)
    subPos = reverse(subPos)

    NewMembership = [sum(masses[ this .⊆ FE ]) for this in subPos]  # Find beliefs
    zPos = range(0,1, length=steps+1)

    zMems = Interval{Float64}[]

    push!(zMems, subPos[1])

    for i = 2:length(zPos)-1                                # Find alpha cuts for return Fuzzy
        theseOnes =  zPos[i] .<= NewMembership .< zPos[i+1]
        if !any(theseOnes);
            thisInt = zMems[end]
        else
            thisInt = hull(subPos[theseOnes])
        end
        push!(zMems, thisInt);
    end

    return Fuzzy(zMems)
end


#=

function DSS2Fuzzy(FE::Array{Interval{T},1}, masses::Array{Float64,1} = ones(Integer(length(FE))) ./ length(FE); steps = length(FE)) where T <: Real

    if length(FE) < steps; steps = length(FE); end

    lefts = PossibilisticArithmetic.left.(FE);
    rights = PossibilisticArithmetic.right.(FE);

    lefts = sort(lefts); rights = sort(rights)   ## Problem ... Would need to mult by 2. But will be over-conservative.

    bins = [lefts, rights];
    bins = unique(bins)

    subPos

    subPos = interval.(lefts, rights)

    subPos = unique(subPos)

    NewMembership = [sum(masses[ this .⊆ FE ]) for this in subPos]  # Find beliefs

    memberships = range(0, 1, length = steps + 1)[1:end - 1]

    these = [findall(α .<= NewMembership) for α in memberships]

    return Fuzzy([hull(subPos[this]) for this in these])
end


function DSS2Fuzzy(FE::Array{Interval{T},1}, masses::Array{Float64,1} = ones(Integer(length(FE))) ./ length(FE); steps = length(FE)) where T <: Real

    if length(FE) < steps; steps = length(FE); end

    lefts = PossibilisticArithmetic.left.(FE);
    rights = PossibilisticArithmetic.right.(FE);

    lefts = sort(lefts); rights = sort(rights, rev = true)   ## Problem ... Would need to mult by 2. But will be over-conservative.

    subPos = interval.(lefts, rights)

    subPos = unique(subPos)

    NewMembership = [sum(masses[ this .⊆ FE ]) for this in subPos]  # Find beliefs

    memberships = range(0, 1, length = steps + 1)[2:end]

    these = [findall(α .<= NewMembership) for α in memberships]

    return Fuzzy([hull(subPos[this]) for this in these])
end
=#

#=
function DSS2Fuzzy(FE::Array{Interval{T},1}, masses::Array{Float64,1} = ones(Integer(length(FE))) ./ length(FE); steps = length(FE)) where T <: Real

    if length(FE) < steps; steps = length(FE); end

    lefts = PossibilisticArithmetic.left.(FE);
    rights = PossibilisticArithmetic.right.(FE);

    lefts = sort(lefts); rights = sort(rights, rev = true)      ## Problem ... Would need to mult by 2. But will be over-conservative.

    subPos = interval.(lefts, rights)

    subPos = unique(subPos)

    NewMembership = [sum(masses[ this .⊆ FE ]) for this in subPos]  # Find beliefs

    zPos = range(0,1, length=steps+1)

    zMems = Interval{Float64}[]
    push!(zMems, subPos[1])

    for i = 2:length(zPos)-1                                # Find alpha cuts for return Fuzzy
        theseOnes =  zPos[i] .<= NewMembership .< zPos[i+1]
        if !any(theseOnes);
            thisInt = zMems[end]
        else
            thisInt = hull(subPos[theseOnes])
        end
        push!(zMems, thisInt);
    end
    return Fuzzy(zMems)
end
=#
function makeCons1(x::Vector{Interval{T}}) where T <: Real

    x = sort(x, lt = ⊂); # x = reverse(x)
    z = Interval{Float64}[]
    for i = 1:length(x) - 1
        thisZ = x[i]
        if !(x[i] ⊆ x[i + 1]); thisZ = hull(x[i], x[i + 1]);end
        push!(z, thisZ)
    end
    push!(z, x[end])
    return reverse(z)
end

dia(x, y) = diam(x) < diam(y)

function makeCons2(x::Vector{Interval{T}}) where T <: Real

    x = sort(x, rev = true);
    z = Interval{Float64}[]
    for i = 1:length(x) - 1
        lo = x[i].lo; hi = x[i].hi;
        if x[i].lo > x[i + 1].lo; lo = x[i + 1].lo;end
        if x[i].hi < x[i + 1].hi; hi = x[i + 1].hi;end
        push!(z, interval(lo, hi))
    end
    push!(z, x[end])
    return z
end

function env(x::FuzzyNumber, y::FuzzyNumber)

    Mems = hull.(x.Membership, y.Membership)
    return FuzzyNumber(Mems)

end


function makepbox(x::FuzzyNumber)

    lefts = left.(x.Membership)
    rights = right.(x.Membership)

    newRights = reverse(rights)

    return pbox(lefts, newRights)
end


function makefuzzy(x::pbox)

    # Need to check if possible

    lefts = x.u;
    rights = x.d;

    newRight = reverse(rights);

    num = Int(floor(length(lefts) / 2));

    Fuzz = FuzzyNumber(interval.(lefts[1:num], newRight[1:num]))

    return descritize(Fuzz, 200)
end

function makefuzzy(x::Union{Real,Interval}, steps = 200)
    Mems = [interval(x) for i = 1:steps]
    return Fuzzy(Mems)
end

function ecdf2fuzzy(x)
    xSort = sort(x)
    l = length(x)
    if iseven(l)
        mid = Integer(l / 2)
        lefts = xSort[1:mid]
        rights = reverse(xSort[mid + 1:end])
        Mems = interval.(lefts, rights)
    else
        mid = Integer((l - 1) / 2)
        lefts = xSort[1:mid]
        rights = reverse(xSort[mid + 2:end])
        core = Interval(xSort[mid + 1])
        Mems = interval.(lefts, rights)
        push!(Mems, core)
    end
    return Fuzzy(Mems)
end

###
#   Checks wether a distribution is inside a fuzzy number
###

function check_inside(F :: Fuzzy, Dist; print = false, usemass = false)
    F_α = F.Membership      # α-cuts
    Numel = length(F_α)
    if print; println("i     |     Nec_α    |     m_dist    "); end
    for (i, el) in  enumerate(F_α)
        m_dist = cdf(Dist, el.hi) - cdf(Dist, el.lo)    # mass of Dist in el
        m_dist = interval(m_dist)
        Nec_α = usemass ? mass(F, el, el).lo : (Numel - i + 1)/Numel
        if print; println("$i     |     $Nec_α    |     $(m_dist.hi)     " ); end
        if Nec_α > m_dist.hi
            return false                                # If nec is higher return false
        end
    end
    return true
end

###
#   Linearly interpolates between Range and Core.
#
#   Downward rounding for lower bound, upward for upper
###

function linearInterp(Core::Interval{T}, Range::Interval{T}, steps = 200)  where T <: Real

    lows = setrounding(T, RoundDown) do
        collect(range(Range.lo, Core.lo, length = steps))
    end

    highs = setrounding(T, RoundUp) do
        collect(range(Range.hi, Core.hi, length = steps))
    end

    return interval.(lows, highs)
end
#=
An alternative to above due to @kolaru

function linearInterp(Core::Interval, Range::Interval, steps = 200)
    width = Range - Core
    return reverse([Core + width * k/steps for k in 1:steps])
end

=#

isfuzzy(x) = typeof(x) <: FuzzyNumber

function Base.show(io::IO, z::FuzzyNumber)

    Range = z.Membership[1]
    Core = z.Membership[end]
    if isscalar(Core); Core = Core.lo; end

    print(io, "Fuzzy: \t ~ ( Range=$Range, Core=$Core )");

end
