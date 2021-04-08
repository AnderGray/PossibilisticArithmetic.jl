

left(x :: IntervalU) = left.(x.v)
right(x :: IntervalU) = right.(x.v)

struct PossNumber <: AbstractPoss

    Core :: IntervalU{T} where T <: Real
    Range :: IntervalU{T} where T <: Real
    Membership :: Vector{IntervalU{T}} where T <: Real

    function PossNumber(Core = intervalU(0.5), Range = intervalU(0, 1), Membership = missing; steps = 200)
        
        if !(Core ⊆ Range); throw(ArgumentError("Core must be a subset of the Range.\nProvided Core = $Core\nProvided Range = $Range")); end

        if ismissing(Membership); 
            Mems = [linearInterp(coreVec, rangeVec, steps) for coreVec in Core.v, rangeVec in Range.v]; 
            reduced = reduce(hcat, Mems)
            Membership = [ ∪(reduced[i,:]) for i =1:steps]
        end

        #if Range != Membership[1]; throw(ArgumentError("Range must be the first element of the Membership function.\nProvided Range = $Range\nProvided Membership[1] = $(Membership[1]))"));end
        #if Core != Membership[end]; throw(ArgumentError("Core must be the last element of the Membership function.\nProvided Core = $Core\nProvided Membership[end] = $(Membership[end]))"));end

        return new(Core, Range, Membership)
    end
end

Possibility = PossNumber

PossNumber(lowerbound :: Real, Core :: Real, upperbound :: Real; steps = 200) = PossNumber(intervalU(Core), intervalU(lowerbound, upperbound), steps = steps)
PossNumber(lowerbound :: Real, CoreLeft :: Real, CoreRight :: Real, upperbound :: Real; steps = 200) = PossNumber(interval(CoreLeft, CoreRight), intervalU(lowerbound, upperbound), steps = steps)
PossNumber(lowerbound :: Real, Core :: Interval{<:Real}, upperbound :: Real; steps = 200) = PossNumber(Core, intervalU(lowerbound, upperbound), steps = steps)

function PossNumber( Membership :: Vector{IntervalU{T}}) where T <: Real

    #if !isnested(Membership); throw(ArgumentError("Invalid Membership function, intervals must be nested from Range (Membership[1]) to Core (Membership[end])")); end
    Core = Membership[end]; Range = Membership[1];
    return PossNumber(Core, Range, Membership, steps = length(Membership))

end

FuzzyNumber(Membership :: Vector{IntervalU{T}}) where T <: Real = PossNumber(Membership)

function membership(F :: PossNumber, x ::Union{Float64, Int64})
    mems = F.Membership
    here = findlast(x .∈ mems)
    if isnothing(here) return 0;end
    return here/(length(mems)+1)
end

function cut( x :: PossNumber, α :: Real )

    if !(α ∈ interval(0, 1) ); throw(ArgumentError("α ∉ [0, 1]\nProvided α = $α")); end

    Mems = deepcopy(x.Membership)
    push!(Mems, Mems[end])

    α = (α *(length( Mems ) -1)) + 1; # Scale α by number of elements
    
    i = Int(round(α, RoundDown))         # Round down for conservatism
    return Mems[i]
end

function Base.show(io::IO, z::PossNumber)

    Range = z.Range
    Core = z.Core
    if length(z.Core.v) == 1; if isscalar(z.Core[1]); Core = z.Core[1].lo; end; end
    
    print(io, "Poss: \t ~ ( Range=$Range, Core=$Core )");

end