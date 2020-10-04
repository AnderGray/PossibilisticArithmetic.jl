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

for op in (:+, :-, :*, :/, :min, :max, :^, :log, :<, :>)
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

function supT(x :: FuzzyNumber, y :: FuzzyNumber; op = +, T = min) 

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


function mobiusTransform2D(x :: Fuzzy, y::Fuzzy, C)

    xMems = x.Membership; yMems = y.Membership;
    xBel  = range(1,0, length=length(xMems))
    yBel  = range(1,0, length=length(yMems))

    cartProd = Array{Interval{Float64},1}[]
    masses = Float64[]

    for i = 1:length(xMems)-1
        for j = 1:length(yMems)-1
            thisMass = C(xBel[i],yBel[j]) - C(xBel[i+1],yBel[j]) - C(xBel[i],yBel[j+1]) + C(xBel[i+1],yBel[j+1])
            push!(masses, thisMass[1])
            push!(cartProd, [xMems[i], yMems[j]])
        end
    end
    return masses, cartProd
end

function sigmaFuzzy(x :: Fuzzy, y::Fuzzy; op = +, C = π())

    if C.func == opp; return levelwiseOpp(x, y, op=op);end
    if C.func == perf; return levelwise(x, y,op= op);end

    zNum = max(length(x.Membership), length(y.Membership))
    masses, cartProd = mobiusTransform2D(x, y, C)
    
    zs = [op(ints[1],ints[2]) for ints in cartProd]

    zUniq = unique(zs)
    
    membership = [sum(masses[ zs .⊇ this ]) for this in zUniq]

    zPos = range(0,1, length=zNum+1)

    zMems = Interval{Float64}[]
    #push!(zMems, zUniq[1])
    
    for i = 1:length(zPos)-1
        theseOnes =  zPos[i] .<= membership .< zPos[i+1]
        if !any(theseOnes); 
            thisInt = zMems[end]
        else
            thisInt = hull(zUniq[theseOnes])
        end
        push!(zMems, thisInt);
    end

    if !isnested(zMems); zMems= makeCons(zMems);end
    return Fuzzy(zMems)
end



function sigmaFuzzy2(x :: Fuzzy, y::Fuzzy; op = +, C = π())

    #if C.func == opp; return levelwiseOpp(x, y, op);end

    zNum = max(length(x.Membership), length(y.Membership))
    masses, cartProd = mobiusTransform2D(x, y, C)
    
    zs = [op(ints[1],ints[2]) for ints in cartProd]

    zUniq = unique(zs)
    
    membership = [sum(masses[ zs .⊇ this ]) for this in zUniq]

    zPos = range(0,1, length=zNum+1)

    p = sortperm(membership)
    membership = membership[p]
    zUniq = zUniq[p];
    zMems = [zUniq[searchsortedfirst(membership, i)] for i in zPos[1:end-1]]
    
    if !isnested(zMems); zMems= makeCons(zMems);end
    return Fuzzy(zMems)
end



function sigmaFuzzy3(x :: Fuzzy, y::Fuzzy; op = +, C = π())

    #if C.func == opp; return levelwiseOpp(x, y, op);end

    zNum = max(length(x.Membership), length(y.Membership))
    masses, cartProd = mobiusTransform2D(x, y, C)
    
    zs = [op(ints[1],ints[2]) for ints in cartProd]

    zUniq = unique(zs)
    
    membership = [sum(masses[ zs .⊇ this ]) for this in zUniq]

    zPos = range(1,0, length=zNum+1)

    p = sortperm(membership)
    membership = membership[p]
    zUniq = zUniq[p];
    zMems = [zUniq[searchsortedlast(membership, i)] for i in zPos[1:end-1]]
    
    push!(zMems, zs[1])

    if !isnested(zMems); zMems= makeCons(zMems);end
    return Fuzzy(zMems)
end



function levelwiseOpp(x::Fuzzy, y :: Fuzzy; op = +)
    xMems = x.Membership; yMems = y.Membership;

    yMems = reverse(yMems)
    zMems = op.(xMems, yMems)

    #zMems = sort(zMems, lt = ⊂, rev = true)
    zMems = makeCons(zMems)
    return Fuzzy(zMems)
end



function tauFuzzy(x :: FuzzyNumber, y :: FuzzyNumber; op = +, C = π())

    xMems = x.Membership; yMems = y.Membership;
    numX = length(xMems); numY = length(yMems);
    xPos  = range(1,0,length=numX); yPos = range(1, 0, length = numY);

    zMems = op.(xMems, yMems);
    zPos  = getindex.(C.(xPos, yPos), 1)

    zPosNew = range(0,1,length=length(zMems)+1)

    zMemsNew = Interval{Float64}[]

    #zRange = zMems[end];
    #push!(zMemsNew, zRange);

    for i = 1:length(zPosNew)-1
        theseOnes =  zPosNew[i] .<= zPos .<= zPosNew[i+1]
        if !any(theseOnes); 
           thisInt = zMemsNew[end]
        else
           thisInt = hull(zMems[theseOnes])
        end
        push!(zMemsNew, thisInt);
    end

    zMemsNew = reverse(zMemsNew)
    return Fuzzy(zMemsNew)


end


function tauFuzzy3(x :: FuzzyNumber, y :: FuzzyNumber; op = +, C = π())

    xMems = x.Membership; yMems = y.Membership;
    numX = length(xMems); numY = length(yMems);
    xBel  = range(1,0,length=numX); yBel = range(1, 0, length = numY);

    zMems = op.(xMems, yMems);
    zPos  = getindex.(C.(xBel, yBel), 1)

    zPosNew = range(0,1,length=length(zMems)+1)

    zMemsNew = Interval{Float64}[]

    zRange = zMems[end];
    #push!(zMemsNew, zRange);

    for i = 1:length(zPosNew)-1
        theseOnes =  zPosNew[i] .<= zPos .<= zPosNew[i+1]
        if !any(theseOnes); 
            thisInt = zMemsNew[end]
        else
            thisInt = hull(zMems[theseOnes])
        end
        push!(zMemsNew, thisInt);
    end
    
    zMemsNew = reverse(zMemsNew)
    #if !isnested(zMemsNew); zMemsNew= makeCons(zMemsNew);end
    return Fuzzy(zMemsNew)

end

function tauFuzzy2(x :: FuzzyNumber, y :: FuzzyNumber; op = +, C = π())

    xMems = x.Membership; yMems = y.Membership;
    numX = length(xMems); numY = length(yMems);
    xBel  = range(1,0,length=numX); yBel = range(1, 0, length = numY);

    zMems = op.(xMems, yMems);
    zPos  = getindex.(C.(xBel, yBel), 1)

    zPosNew = range(0,1,length=length(zMems)+1)

    zMemsNew = Interval{Float64}[]

    zRange = zMems[end];
    #push!(zMemsNew, zRange);

    p = sortperm(zPos)
    zPos = zPos[p]
    zMems = zMems[p];
    zMemsNew = [zMems[searchsortedfirst(zPos, i)] for i in zPosNew[1:end-1]]

    zMemsNew = reverse(zMemsNew)
    if !isnested(zMemsNew); zMemsNew= makeCons(zMemsNew);end
    return Fuzzy(zMemsNew)

end


Tgen(x,y) = min(min(1,2*x), min(1,2*y))
Tind(x,y) = min(1 - (1-x)^2, 1 - (1-y)^2)

TestNorm(x, y) = 2 .* perf(x,y) .- opp(x,y)

supMin(x:: FuzzyNumber, y :: FuzzyNumber;op=+) = supT(x, y, op=op, T = min)
supGen(x:: FuzzyNumber, y :: FuzzyNumber;op=+) = supT(x, y, op=op, T = Tgen)
supInd(x:: FuzzyNumber, y :: FuzzyNumber;op=+) = supT(x, y, op=op, T = Tind)




function supCop(x :: FuzzyNumber, y :: FuzzyNumber; op = +, C = min) 

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
    zPos = [C(x,y)[1] for x in xPos, y in yPos]

    zMem = [zRange for i = 1:zNumMem];  # Begin z membership as just a vector of z's range

    zCoreVal = zs[zPos .== 1]           # Find all values in core
    zCore = interval(minimum(zCoreVal), maximum(zCoreVal));     # Construct core

    zMem[end] = zCore
    zPs = range(0, 1, length = zNumMem)         # alpha values for z

    for i = 2:(zNumMem-1)                       # Iterate through alpha values, and contruct intervals
        zVals = zs[ .!(zPs[i-1] .<= zPos .<= zPs[i+1])]
        zInt = interval(minimum(zVals), maximum(zVals))
        zMem[i] = zInt
    end

    return FuzzyNumber(zMem)

end



function supCop2(x :: FuzzyNumber, y :: FuzzyNumber; op = +, C = min) 

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
    zPos = [C(x,y)[1] for x in xPos, y in yPos]

    #zMem = [zRange for i = 1:zNumMem];  # Begin z membership as just a vector of z's range

    #zCoreVal = zs[zPos .== 1]           # Find all values in core
    #zCore = interval(minimum(zCoreVal), maximum(zCoreVal));     # Construct core

    zMems = range(zRange.lo, zRange.hi, length = zNumMem)


    zPs = zeros(zNumMem)

    for i = 2:zNumMem-1                       # Iterate through alpha values, and contruct intervals
        zPs[i] = maximum(zPos[.!(zMems[i-1] .<= zs .<= zMems[i+1])])
        #zInt = interval(minimum(zVals), maximum(zVals))
        #zMem[i] = zInt
    end

    return zPs, zMems

end