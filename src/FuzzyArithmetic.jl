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

for op in (:-, :sin, :cos, :tan, :exp, :log)
    @eval ($op)( x::FuzzyNumber) = FuzzyNumber(broadcast($op, x.Membership))
end

#=
+( x :: FuzzyNumber, y :: FuzzyNumber)  = levelwise(x, y, op = +)
-( x :: FuzzyNumber, y :: FuzzyNumber)  = levelwise(x, y, op = -)
*( x :: FuzzyNumber, y :: FuzzyNumber)  = levelwise(x, y, op = *)
/( x :: FuzzyNumber, y :: FuzzyNumber)  = levelwise(x, y, op = /)
=#


function supMin(x :: FuzzyNumber, y :: FuzzyNumber; op = +)

    zRange = op(x.Range, y.Range);

    xNumMem = length(x.Membership); yNumMem = length(y.Membership);

    zNumMem = max(xNumMem, yNumMem);    # Number of membeship elements Z

    xLeft = left.(x.Membership); xRight = right.(x.Membership);
    yLeft = left.(y.Membership); yRight = right.(y.Membership);

    xs = [xLeft; reverse(xRight)]; 
    ys = [yLeft; reverse(yRight)];

    xPos = range(0, 1, length = xNumMem);
    yPos = range(0, 1, length = yNumMem);

    xPos = [xPos; reverse(xPos)]
    yPos = [yPos; reverse(yPos)]

    zs   = [map(op, x, y) for x in xs, y in ys]
    zPos = [min(x, y) for x in xPos, y in yPos]

    zMem = [zRange for i = 1:zNumMem];  # Begin z membership as just a vector of z's range

    zCoreVal = zs[zPos .== 1]           # Find all values in core
    zCore = interval(minimum(zCoreVal), maximum(zCoreVal));     # Construct core

    zMem[end] = zCore

    zPs = range(0, 1, length = zNumMem)         # alpha values for z

    for i = 2:(zNumMem-1)                       # Iterate through alpha values, and contruct intervals
        zVals = zs[zPs[i-1] .<= zPos .<= zPs[i+1]]
        zInt = interval(minimum(zVals), maximum(zVals))
        zMem[i] = zInt
    end

    return FuzzyNumber(zMem)

end

function supInd(x :: FuzzyNumber, y :: FuzzyNumber; op = +)

    zRange = op(x.Range, y.Range);

    xNumMem = length(x.Membership); yNumMem = length(y.Membership);

    zNumMem = max(xNumMem, yNumMem);    # Number of membeship elements Z

    xLeft = left.(x.Membership); xRight = right.(x.Membership);
    yLeft = left.(y.Membership); yRight = right.(y.Membership);

    xs = [xLeft; reverse(xRight)]; 
    ys = [yLeft; reverse(yRight)];

    xPos = range(0, 1, length = xNumMem);
    yPos = range(0, 1, length = yNumMem);

    xPos = [xPos; reverse(xPos)]
    yPos = [yPos; reverse(yPos)]

    zs = [map(op, x, y) for x in xs, y in ys]
    zPos = [min(1-(1-x)^2,1-(1-y)^2) for x in xPos, y in yPos]

    zMem = [zRange for i = 1:zNumMem];  # Begin z membership as just a vector of z's range

    zCoreVal = zs[zPos .== 1]           # Find all values in core
    zCore = interval(minimum(zCoreVal), maximum(zCoreVal));     # Construct core

    zMem[end] = zCore

    zPs = range(0, 1, length = zNumMem)         # alpha values for z

    for i = 2:(zNumMem-1)                       # Iterate through alpha values, and contruct intervals
        zVals = zs[zPs[i-1] .<= zPos .<= zPs[i+1]]
        zInt = interval(minimum(zVals), maximum(zVals))
        zMem[i] = zInt
    end

    return FuzzyNumber(zMem)

    
end

#Tgen(x,y) = min(min(1,2*x), min(1,2*y))
#Tindp(x,y) = min(1 - (1-x)^2, 1 - (1-y)^2)

function supGen(x :: FuzzyNumber, y :: FuzzyNumber; op = +)

    zRange = op(x.Range, y.Range);

    xNumMem = length(x.Membership); yNumMem = length(y.Membership);

    zNumMem = max(xNumMem, yNumMem);    # Number of membeship elements Z

    xLeft = left.(x.Membership); xRight = right.(x.Membership);
    yLeft = left.(y.Membership); yRight = right.(y.Membership);

    xs = [xLeft; reverse(xRight)]; 
    ys = [yLeft; reverse(yRight)];

    xPos = range(0, 1, length = xNumMem);
    yPos = range(0, 1, length = yNumMem);

    xPos = [xPos; reverse(xPos)]
    yPos = [yPos; reverse(yPos)]

    zs = [map(op, x, y) for x in xs, y in ys]
    zPos = [min(min(1, 2*x), min(1, 2*y)) for x in xPos, y in yPos]

    zMem = [zRange for i = 1:zNumMem];  # Begin z membership as just a vector of z's range

    zCoreVal = zs[zPos .== 1]           # Find all values in core
    zCore = interval(minimum(zCoreVal), maximum(zCoreVal));     # Construct core

    zMem[end] = zCore

    zPs = range(0, 1, length = zNumMem)         # alpha values for z

    for i = 2:(zNumMem-1)                       # Iterate through alpha values, and contruct intervals
        zVals = zs[zPs[i-1] .<= zPos .<= zPs[i+1]]
        zInt = interval(minimum(zVals), maximum(zVals))
        zMem[i] = zInt
    end

    return FuzzyNumber(zMem)


end

function supT(x :: FuzzyNumber, y :: FuzzyNumber; op = +, T = M():: tnorm) 

    zRange = op(x.Range, y.Range);

    xNumMem = length(x.Membership); yNumMem = length(y.Membership);

    zNumMem = max(xNumMem, yNumMem);    # Number of membeship elements Z

    xLeft = left.(x.Membership); xRight = right.(x.Membership);
    yLeft = left.(y.Membership); yRight = right.(y.Membership);

    xs = [xLeft; reverse(xRight)]; 
    ys = [yLeft; reverse(yRight)];

    xPos = range(0, 1, length = xNumMem);
    yPos = range(0, 1, length = yNumMem);

    xPos = [xPos; reverse(xPos)]
    yPos = [yPos; reverse(yPos)]

    zs = [map(op, x, y) for x in xs, y in ys]
    zPos = [T(x,y)[1] for x in xPos, y in yPos]

    zMem = [zRange for i = 1:zNumMem];  # Begin z membership as just a vector of z's range

    zCoreVal = zs[zPos .== 1]           # Find all values in core
    zCore = interval(minimum(zCoreVal), maximum(zCoreVal));     # Construct core

    zMem[end] = zCore

    zPs = range(0, 1, length = zNumMem)         # alpha values for z

    for i = 2:(zNumMem-1)                       # Iterate through alpha values, and contruct intervals
        zVals = zs[zPs[i-1] .<= zPos .<= zPs[i+1]]
        zInt = interval(minimum(zVals), maximum(zVals))
        zMem[i] = zInt
    end

    return FuzzyNumber(zMem)

end