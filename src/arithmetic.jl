###
#   This file is a part of PossibilisticArithmetic.jl package
#
#   Defines fuzzy arithemtic functions
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
###


# For now we will use abstract type. When library is finished will use concrete
function levelwise(x::AbstractPoss, y::AbstractPoss; op = +)

    numX = length(x.Membership); numY = length(y.Membership);

    if numX != numY
        num = max(numX, numY)
        x = descritize(x, num); y = descritize(y, num)
    end

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

function levelwiseOpp(x::FuzzyNumber, y::FuzzyNumber; op = +)

    numX = length(x.Membership); numY = length(y.Membership);

    if numX != numY
        num = max(numX, numY)
        x = descritize(x, num); y = descritize(y, num)
    end

    xMems = x.Membership; yMems = y.Membership;

    yMems = reverse(yMems)
    zMems = op.(xMems, yMems)

    return DSS2Fuzzy(zMems)
end


for op in (:+, :-, :*, :/, :min, :max, :^, :log, :<, :>)
    @eval ($op)(x::AbstractPoss, y::AbstractPoss) = levelwise(x, y, op = $op)
    @eval ($op)(x::AbstractPoss, n::Real) = FuzzyNumber(broadcast($op, x.Membership, n))
    @eval ($op)(n::Real, x::AbstractPoss) = FuzzyNumber(broadcast($op, n, x.Membership))
end

for op in (:-, :sin, :cos, :tan, :exp, :log)
    @eval ($op)(x::AbstractPoss) = FuzzyNumber(broadcast($op, x.Membership))
end

#=
+( x :: FuzzyNumber, y :: FuzzyNumber)  = levelwise(x, y, op = +)
-( x :: FuzzyNumber, y :: FuzzyNumber)  = levelwise(x, y, op = -)
*( x :: FuzzyNumber, y :: FuzzyNumber)  = levelwise(x, y, op = *)
/( x :: FuzzyNumber, y :: FuzzyNumber)  = levelwise(x, y, op = /) =#

function mobiusTransform2D(x::Fuzzy, y::Fuzzy, C)

    xMems = x.Membership; yMems = y.Membership;
    xBel  = range(1, 0, length = length(xMems)+1)
    yBel  = range(1, 0, length = length(yMems)+1)

    cartProd = Array{Interval{Float64},1}[]
    masses = Float64[]

    for i = 1:length(xMems)
        for j = 1:length(yMems)
            thisMass = C(xBel[i], yBel[j]) - C(xBel[i + 1], yBel[j]) - C(xBel[i], yBel[j + 1]) + C(xBel[i + 1], yBel[j + 1])
            push!(masses, thisMass[1])
            push!(cartProd, [xMems[i], yMems[j]])
        end
    end
    return masses, cartProd
end

function sigmaFuzzy(x::Fuzzy, y::Fuzzy; op = +, C = Pi())


    if op == -;
        return sigmaFuzzy(x, -y, op=+, C = rotate270(C));
    end
    if op == /;
        return sigmaFuzzy(x, 1/y, op=/, C = rotate270(C));
    end

    if C == M
        if op == -; return levelwiseOpp(x, -y, op = +); end
        if op== / ; return levelwiseOpp(x, 1/y, op = +); end
        return levelwise(x, y, op = op)
    end

    if C == W
        if op == -; return levelwise(x, -y, op = +); end
        if op== / ; return levelwise(x, 1/y, op = +);end
        return levelwiseOpp(x, y, op = op)
    end


    zNum = max(length(x.Membership), length(y.Membership))
    masses, cartProd = mobiusTransform2D(x, y, C)       # Get carteesian prod and masses from Möbius

    zs = [op(ints[1], ints[2]) for ints in cartProd]     # Evaluate carteesian product with interval arithm

    return DSS2Fuzzy(zs, masses, steps = max(length(x.Membership), length(y.Membership)))
end

sigmaFuzzy(x::FuzzyNumber, y::Real; op = +, C = π()) = sigmaFuzzy(x, makefuzzy(y); op = op, C = C)
sigmaFuzzy(x::Real, y::FuzzyNumber; op = +, C = π()) = sigmaFuzzy(makefuzzy(x), y; op = op, C = C)

function sigmaFuzzySampling(x::Fuzzy, y::Fuzzy; op = +, C = Pi(), Nsamps = 10^3)

    if op == -;
        return sigmaFuzzySampling(x, -y, op=+, C = rotate270(C), Nsamps);
    end
    if op == /;
        return sigmaFuzzySampling(x, 1/y, op=/, C = rotate270(C), Nsamps);
    end

    masses, cartProd = mobiusTransform2D(x, y, C)       # Get carteesian prod and masses from Möbius

    x1 = x.Membership[1]; y1 = y.Membership[1];

    xSamps = rand(Nsamps) .* diam(x1) .+ x1.lo
    ySamps = rand(Nsamps) .* diam(y1) .+ y1.lo

    #XYsamps = hcat(xSamps, ySamps);
    XYsamps = IntervalBox.(xSamps, ySamps);
    boxes = IntervalBox.(cartProd);

    outSamps = op(xSamps, ySamps);

    outInts = interval.(zeros(length(masses)));
    these = [samp ⊆  b for b in boxes, samp in XYsamps];

    those = []

    for i =1:length(masses)
        if !any(these[i,:])
            push!(those, i)
        else
            outInts[i] = hull(interval.(outSamps[these[i,:]]))
        end
    end

    those = sort!(those, rev= true)

    [popat!(masses, t) for t in those]
    [popat!(outInts, t) for t in those]

    return DSS2Fuzzy(outInts, masses, steps = max(length(x.Membership), length(y.Membership)))
end

function sigmaFuzzySampling2(x::Fuzzy, y::Fuzzy; op = +, C = Pi(), Nsamps =10^3)


    masses, cartProd = mobiusTransform2D(x, y, C)       # Get carteesian prod and masses from Möbius

    x1 = x.Membership[1]; y1 = y.Membership[1];

    xSamps = rand(Nsamps) .* diam(x1) .+ x1.lo
    ySamps = rand(Nsamps) .* diam(y1) .+ y1.lo

    #XYsamps = hcat(xSamps, ySamps);
    XYsamps = IntervalBox.(xSamps, ySamps);
    boxes = IntervalBox.(cartProd);

    outSamps = op(xSamps, ySamps);

    poss = zeros(Nsamps)

    for i = 1:Nsamps
        these = XYsamps[i,:] .⊆ boxes
        poss[i] = sum(masses[these])
    end

    #poss = 1 .- poss
    return outSamps, poss
end


function sigmaFuzzySampling3(x::Fuzzy, y::Fuzzy; op = +, C = Pi(), Nsamps =10^3)


    masses, cartProd = mobiusTransform2D(x, y, C)       # Get carteesian prod and masses from Möbius

    x1 = x.Membership[1]; y1 = y.Membership[1];

    xSamps = rand(Nsamps) .* diam(x1) .+ x1.lo
    ySamps = rand(Nsamps) .* diam(y1) .+ y1.lo

    XYsamps = hcat(xSamps, ySamps);
    #XYsamps = IntervalBox.(xSamps, ySamps);
    boxes = IntervalBox.(cartProd);

    outSamps = op(xSamps, ySamps);

    poss = zeros(Nsamps)

    for i = 1:Nsamps
        these = [XYsamps[i,:] ∈ b for b in boxes]
        poss[i] = sum(masses[these])
    end

    #poss = 1 .- poss
    return outSamps, poss
end

function sigmaFuzzySampling4(x::Fuzzy, y::Fuzzy; op = +, C = Pi(), Nsamps =10^3)


    masses, cartProd = mobiusTransform2D(x, y, C)       # Get carteesian prod and masses from Möbius

    x1 = x.Membership[1]; y1 = y.Membership[1];

    xSamps = rand(Nsamps) .* diam(x1) .+ x1.lo
    ySamps = rand(Nsamps) .* diam(y1) .+ y1.lo

    #XYsamps = hcat(xSamps, ySamps);
    XYsamps = IntervalBox.(xSamps, ySamps);
    boxes = IntervalBox.(cartProd);

    subposs = x.Membership .× y.Membership

    #=
    verts = vertices.(boxes)
    ver = verts[1]
    for i = 2:length(boxes)-1
        ver = vcat(ver, verts[i])
    end
    =#

    plaus = zeros(length(subposs))
    for i = 1:length(subposs)
        these = [subposs[i] ⊆ b for b in boxes]
        plaus[i] = sum(masses[these])
    end

    outSamps = op(xSamps, ySamps);

    poss = zeros(Nsamps)

    for i = 1:Nsamps
        these = XYsamps[i,:] .⊆ subposs
        poss[i] = maximum(plaus[these])
    end


    return outSamps, poss
end

function vertices(box :: IntervalBox)

    Ndims = length(box);
    Nverts = 2^Ndims;

    verts = zeros(Nverts, Ndims);

    for i =1:Ndims

        A = [box[i].lo box[i].hi]
        j1 = 2^(i-1); j2 = 2^(Ndims - i);

        As = repeat(A , outer = (j1, j2));

        verts[:,i] = As[:]

    end

    return verts

end

function tauFuzzy(x::FuzzyNumber, y::FuzzyNumber; op = +, C = Pi())

    if op == -; return tauFuzzy(x, -y, op=+, C = rotate270(C)); end
    if op == /; return tauFuzzy(x, 1/y, op=/, C = rotate270(C)); end

    numX = length(x.Membership); numY = length(y.Membership);

    if numX != numY
        num = max(numX, numY)
        x = descritize(x, num); y = descritize(y, num)
    end

    xMems = x.Membership; yMems = y.Membership;
    numX = length(xMems); numY = length(yMems);
    xBel  = range(1, 0, length = numX); yBel = range(1, 0, length = numY);

    zMems = op.(xMems, yMems);
    zBel  = getindex.(C.(xBel, yBel), 1)

    zPosNew = range(0, 1, length = length(zMems) + 1)

    zMemsNew = Interval{Float64}[]

    for i = 1:length(zPosNew) - 1
        theseOnes =  zPosNew[i] .<= zBel .<= zPosNew[i + 1]
        if !any(theseOnes);
            thisInt = zMemsNew[end]
        else
            thisInt = hull(zMems[theseOnes])
        end
        push!(zMemsNew, thisInt);
    end

    zMemsNew = reverse(zMemsNew)
    return FuzzyNumber(zMemsNew)


end

tauFuzzy(x::FuzzyNumber, y::Real; op = +, C = π()) = tauFuzzy(x, makefuzzy(y); op = +, C = Pi())
tauFuzzy(x::Real, y::FuzzyNumber; op = +, C = π()) = tauFuzzy(makefuzzy(x), y; op = +, C = Pi())



##
#   OldSchool Williamson & Downs T-norm convolution
##
function supT(x::FuzzyNumber, y::FuzzyNumber; op = +, T = min)

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
    zPos = [T(x, y)[1] for x in xPos, y in yPos]

    zMem = [zRange for i = 1:zNumMem];  # Begin z membership as just a vector of z's range

    zCoreVal = zs[zPos .== 1]           # Find all values in core
    zCore = interval(minimum(zCoreVal), maximum(zCoreVal));     # Construct core

    zMem[end] = zCore

    zPs = range(0, 1, length = zNumMem)         # alpha values for z

    for i = 2:(zNumMem - 1)                       # Iterate through alpha values, and contruct intervals
        zVals = zs[zPs[i - 1] .<= zPos .<= zPs[i + 1]]
        zInt = interval(minimum(zVals), maximum(zVals))
        zMem[i] = zInt
    end

    return FuzzyNumber(zMem)

end

Tgen(x,y) = min(min(1, 2 * x), min(1, 2 * y))
Tind(x,y) = min(1 - (1 - x)^2, 1 - (1 - y)^2)

TestNorm(x, y) = 2 .* perf(x, y) .- opp(x, y)

supMin(x::FuzzyNumber, y::FuzzyNumber;op = +) = supT(x, y, op = op, T = min)
supGen(x::FuzzyNumber, y::FuzzyNumber;op = +) = supT(x, y, op = op, T = Tgen)
supInd(x::FuzzyNumber, y::FuzzyNumber;op = +) = supT(x, y, op = op, T = Tind)




function sigmaFuzzyOld(x :: Fuzzy, y::Fuzzy; op = +, C = Pi())

    if C == W; return levelwiseOpp(x, y, op=op);end
    if C == M; return levelwise(x, y, op= op);end

    zNum = max(length(x.Membership), length(y.Membership))
    masses, cartProd = mobiusTransform2D(x, y, C)       # Get carteesian prod and masses from Möbius

    zs = [op(ints[1],ints[2]) for ints in cartProd]     # Evaluate carteesian product with interval arithm

    zUniq = unique(zs)                                  # Remove repeated intervals for efficiency

    #=
    lefts = left.(zUniq);
    rights = sort(right.(zUniq),rev = true);

    subPos = interval.(lefts, rights)
    =#

    membership = [sum(masses[ zs .⊇ this ]) for this in zUniq]  # Find beliefs

    zPos = range(0,1, length=zNum+1)

    zMems = Interval{Float64}[]
    push!(zMems, zUniq[1])

    for i = 2:length(zPos)-1                                # Find alpha cuts for return Fuzzy
        theseOnes =  zPos[i] .<= membership .< zPos[i+1]
        if !any(theseOnes);
            thisInt = zMems[end]
        else
            thisInt = hull(zUniq[theseOnes])
        end
        push!(zMems, thisInt);
    end
    #zMems = sort(zMems, lt = ⊂, rev= true)
    if !isnested(zMems); zMems= makeCons(zMems);end
    return Fuzzy(zMems)
end
#=  Old codess

function sigmaFuzzy1(x :: Fuzzy, y::Fuzzy; op = +, C = Pi())

    #if C.func == opp; return levelwiseOpp(x, y, op);end

    zNum = max(length(x.Membership), length(y.Membership))
    masses, cartProd = mobiusTransform2D(x, y, C)       # Get carteesian prod and masses from Möbius

    zs = [op(ints[1],ints[2]) for ints in cartProd]     # Evaluate carteesian product with interval arithm

    zUniq = unique(zs)                                  # Remove repeated intervals for efficiency

    lefts = left.(zUniq);
    rights = sort(right.(zUniq), rev= true);

    subPos = interval.(lefts, rights)

    membership = [sum(masses[ zs .⊇ this ]) for this in subPos]  # Find beliefs

    zPos = range(0,1, length=zNum+1)

    p = sortperm(membership)                                    # sort according to membership
    membership = membership[p]
    zUniq = zUniq[p];
    zMems = [zUniq[searchsortedfirst(membership, i)] for i in zPos[1:end-1]]    # Find alpha cuts for return Fuzzy

    if !isnested(zMems); zMems= makeCons(zMems);end
    return Fuzzy(zMems)
end


function sigmaFuzzy3(x :: FuzzyNumber, y::FuzzyNumber; op = +, C = Pi())

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




function tauFuzzy3(x :: FuzzyNumber, y :: FuzzyNumber; op = +, C = Pi())

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

function tauFuzzy2(x :: FuzzyNumber, y :: FuzzyNumber; op = +, C = Pi())

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



function supCop(x :: FuzzyNumber, y :: FuzzyNumber; op = +, C = min)

    zRange = op(x.Range, y.Range);

    xNumMem = length(x.Membership); yNumMem = length(y.Membership);

    zNumMem = max(xNumMem, yNumMem);    # Number of membeship elements Z

    xLeft = left.(x.Membership); xRight = right.(x.Membership);
    yLeft = left.(y.Membership); yRight = right.(y.Membership);

    xs = [xLeft; reverse(xRight)];
    ys = [yLeft; reverse(yRight)];

    xPos = range(1, 0, length = xNumMem);
    yPos = range(1, 0, length = yNumMem);

    xPos = [xPos; reverse(xPos)]
    yPos = [yPos; reverse(yPos)]

    zs = [map(op, x, y) for x in xs, y in ys]
    zPos = [1 - C(x,y)[1] for x in xPos, y in yPos]

    zMem = [zRange for i = 1:zNumMem];  # Begin z membership as just a vector of z's range

    zCoreVal = zs[zPos .== 1]           # Find all values in core
    zCore = interval(minimum(zCoreVal), maximum(zCoreVal));     # Construct core

    zMem[end] = zCore
    zPs = range(0, 1, length = zNumMem)         # alpha values for z

    for i = 2:(zNumMem-1)                       # Iterate through alpha values, and contruct intervals
        zVals = zs[(zPs[i-1] .<= zPos .<= zPs[i+1])]
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

end =#
