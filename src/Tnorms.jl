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

abstract type AbstractTnorm <: Real end

mutable struct tnorm <: AbstractTnorm

    T ::  Array{Float64,2}
    func :: Union{Function, Missing}
    param :: Union{Float64, Missing}

    function tnorm(T = missing; func = missing, param = missing)

        #if (ismissing(cdf) && ismissing(density)) throw(ArgumentError("Cdf and/or density must be provided")) end
        if (ismissing(T));throw(ArgumentError("Cdf and/or density must be provided")); end
        return new(T, func, param)
    end
end

function (obj :: tnorm)(X, Y, useInterp = false)             # Evaluating the copula

    #if (!issorted(X) || !issorted(Y)) throw(ArgumentError("cdf evaluation request must be given in assending order"));end

    if !useInterp
        if !ismissing(obj.func)
            if !ismissing(obj.param)
                return obj.func(X,Y,obj.param)
            end
            return obj.func(X,Y)
        end
    end

    xsize, ysize = size(obj.T)

    xIndexLower = Int.(floor.(X .* (xsize-1)) .+ 1)
    yIndexLower = Int.(floor.(Y .* (ysize-1)) .+ 1)

    xIndexUpper = Int.(ceil.(X .* (xsize-1)) .+ 1)
    yIndexUpper = Int.(ceil.(Y .* (ysize-1)) .+ 1)

    lowerVals = getindex(obj.T, xIndexLower, yIndexLower)
    upperVals = getindex(obj.T, xIndexUpper, yIndexUpper)
    return interval.(lowerVals, upperVals);
end



indep(X, Y) = [x*y for x in X, y in Y];
perf(X, Y)  = [min(x,y) for x in X, y in Y];
opp(X, Y)   = [max(x+y-1,0) for x in X, y in Y];
zee(X,Y)    = [x == 1.0 ? y : y == 1.0 ? x : 0.0  for x in X, y in Y ]

coindep(X, Y) = 1 .- [(1-x)*(1-y) for x in X, y in Y];
coperf(X, Y)  = 1 .- [min((1-x),(1-y)) for x in X, y in Y];
coopp(X, Y)   = 1 .- [max((1-x)+(1-y)-1,0) for x in X, y in Y];
cozee(X,Y)    = 1 .- [(1-x) == 1.0 ? y : (1-x) == 1.0 ? x : 0.0  for x in X, y in Y ]

function π( steps = 200 )
    x = range(0, 1, length = steps)
    return tnorm(indep(x, x), func = indep)
end

function M( steps = 200 )
    x = range(0, 1, length = steps)
    return tnorm(perf(x, x), func = perf)
end

function W( steps = 200 )
    x = range(0, 1, length = steps)
    return tnorm(opp(x, x), func = opp)
end

function Z( steps = 200 )
    x = range(0, 1, length = steps)
    return tnorm(zee(x, x), func = zee)
end

function coπ( steps = 200 )
    x = range(0, 1, length = steps)
    return tnorm(coindep(x, x), func = coindep)
end

function coM( steps = 200 )
    x = range(0, 1, length = steps)
    return tnorm(coperf(x, x), func = coperf)
end

function coW( steps = 200 )
    x = range(0, 1, length = steps)
    return tnorm(coopp(x, x), func = coopp)
end

function coZ( steps = 200 )
    x = range(0, 1, length = steps)
    return tnorm(cozee(x, x), func = cozee)
end

function coZ( steps = 200 )
    x = range(0, 1, length = steps)
    return tnorm(cozee(x, x), func = cozee)
end


function Base.show(io::IO, z::tnorm)

    statement1 = "Arbitrary"
    statement2 = ""

    if (!ismissing(z.func))
        func = z.func
        if (func == indep); func ="π";end
        if (func == perf); func = "M";end
        if (func == opp); func = "W";end
        if (func == zee); func = "Z";end
        if (func == coindep); func ="coπ";end
        if (func == coperf); func = "coM";end
        if (func == coopp); func = "coW";end
        if (func == cozee); func = "coZ";end
        statement1 = "$(func)"
    end

    if (!ismissing(z.param))
        func = z.func
        parName = "par"
        statement2 = "$parName=$(z.param)"
    end

    print(io, "T-norm ~ $statement1($statement2)");
end