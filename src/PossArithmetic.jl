


function levelwise(x::AbstractPoss, y::AbstractPoss; op = +)

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

    return PossNumber(memZ)

end

function levelwiseOpp(x::PossNumber, y::PossNumber; op = +)
    xMems = x.Membership; yMems = y.Membership;

    yMems = reverse(yMems)
    zMems = op.(xMems, yMems)

    zMems = reverse(zMems)
    return PossNumber(zMems)
end




function mobiusTransform2D(x::PossNumber, y::PossNumber, C)

    xMems = x.Membership; yMems = y.Membership;
    xBel  = range(1, 0, length = length(xMems))
    yBel  = range(1, 0, length = length(yMems))

    cartProd = Array{IntervalU{Float64},1}[]
    masses = Float64[]

    for i = 1:length(xMems) - 1
        for j = 1:length(yMems) - 1
            thisMass = C(xBel[i], yBel[j]) - C(xBel[i + 1], yBel[j]) - C(xBel[i], yBel[j + 1]) + C(xBel[i + 1], yBel[j + 1])
            push!(masses, thisMass[1])
            push!(cartProd, [xMems[i], yMems[j]])
        end
    end
    return masses, cartProd
end

function sigmaFuzzy(x::PossNumber, y::PossNumber; op = +, C = Pi())

    if C.func == BivariateCopulas.opp; return levelwiseOpp(x, y, op = op);end
    if C.func == BivariateCopulas.perf; return levelwise(x, y, op = op);end

    zNum = max(length(x.Membership), length(y.Membership))
    masses, cartProd = mobiusTransform2D(x, y, C)       # Get carteesian prod and masses from Möbius

    zs = [op(ints[1], ints[2]) for ints in cartProd]     # Evaluate carteesian product with interval arithm

    # zUniq = unique(zs)                                  # Remove repeated intervals for efficiency
    zUniq = zs
    #= 
    lefts = left.(zUniq);
    rights = sort(right.(zUniq),rev = true);

    subPos = interval.(lefts, rights) =#

    mems = [sum(masses[ zs .⊇ this ]) for this in zUniq]  # Find beliefs

    zPos = range(0, 1, length = zNum + 1)

    zMems = IntervalU{Float64}[]
    push!(zMems, zUniq[1])

    for i = 2:length(zPos) - 1                                # Find alpha cuts for return Fuzzy
        theseOnes =  zPos[i] .<= mems .< zPos[i + 1]
        if !any(theseOnes);
            thisInt = zMems[end]
        else
            sets = reduce(vcat, getfield.(zs[theseOnes], :v))
            thisInt = ∪(sets)
        end
        push!(zMems, thisInt);
    end
    # zMems = sort(zMems, lt = ⊂, rev= true)
    # if !isnested(zMems); zMems= makeCons(zMems);end
    return PossNumber(zMems)
end


function tauFuzzy(x::AbstractPoss, y::AbstractPoss; op = +, C = Pi())

    xMems = x.Membership; yMems = y.Membership;
    numX = length(xMems); numY = length(yMems);
    xBel  = range(1, 0, length = numX); yBel = range(1, 0, length = numY);

    zMems = op.(xMems, yMems);
    zBel  = getindex.(C.(xBel, yBel), 1)

    zPosNew = range(0, 1, length = length(zMems) + 1)

    zMemsNew = IntervalU{Float64}[]

    # zRange = zMems[end];
    # push!(zMemsNew, zRange);

    for i = 1:length(zPosNew) - 1
        theseOnes =  zPosNew[i] .<= zBel .<= zPosNew[i + 1]
        if !any(theseOnes);
            thisInt = zMemsNew[end]
        else
            sets = reduce(vcat, getfield.(zMems[theseOnes], :v))
            thisInt = ∪(sets)
        end
        push!(zMemsNew, thisInt);
    end

    zMemsNew = reverse(zMemsNew)
    return PossNumber(zMemsNew)
end


∩(x::PossNumber, y::PossNumber) = PossNumber(.∩(x.Membership, y.Membership)[1:100])
