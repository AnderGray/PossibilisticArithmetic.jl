###
#   This file is a part of FuzzyArithmetic.jl package
#
#   Defines fuzzy arithemtic functions
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
###


# For now we will use abstract type. When library is finished will use concrete
function levelwise(x :: FuzzyNumber, y :: FuzzyNumber; op = +)

    memX = x.Membership; memY = y.Membership;

    lenX = length(memX); lenY = length(memY);

    if lenX == lenY
        memZ = map(op, memX, memY)
    elseif lenX > lenY
        memY = cut.(y, range(0, 1, length = lenX))
        memZ = map(op, memX, memY)
    elseif lenX < lenY
        memX = cut.(x, range(0, 1, length = lenY))
        memZ = map(op, memX, memY)
    end

    return FuzzyNumber(memZ)

end

##
#   https://github.com/JuliaStats/StatsBase.jl/blob/6f9952d5c9a92faa851822b2e06051f7fc59c294/src/hist.jl#L226-L230
##
for op in (:+, :-, :*, :/, :min, :max, :^, :log)
    @eval ($op)( x::FuzzyNumber, y::FuzzyNumber) = levelwise(x, y, op = $op)
    @eval ($op)( x::FuzzyNumber, n::Real) = FuzzyNumber(broadcast($op, x.Membership, n))
    @eval ($op)( n::Real, x::FuzzyNumber) = FuzzyNumber(broadcast($op, n, x.Membership))
end

for op in (:sin, :cos, :tan, :exp, :log)
    @eval ($op)( x::FuzzyNumber) = FuzzyNumber(broadcast($op, x.Membership))
end

#=
+( x :: FuzzyNumber, y :: FuzzyNumber)  = levelwise(x, y, op = +)
-( x :: FuzzyNumber, y :: FuzzyNumber)  = levelwise(x, y, op = -)
*( x :: FuzzyNumber, y :: FuzzyNumber)  = levelwise(x, y, op = *)
/( x :: FuzzyNumber, y :: FuzzyNumber)  = levelwise(x, y, op = /)
=#


#function maxMin(x :: FuzzyNumber, y :: FuzzyNumber; op = +); end

#function supT(x :: FuzzyNumber, y :: FuzzyNumber; op = +); end