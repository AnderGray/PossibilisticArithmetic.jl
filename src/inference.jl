###
#   This file is a part of FuzzyArithmetic.jl package
#
#   Functions for performing inverse problems with fuzzy numbers
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
###


# Note, not riggerous
function invertUniv(f :: Function, outFuzzys :: FuzzyNumber, InputRange, Nsamples = 10^4, Ndes = 200)

    samps = rand(Nsamples) .* (InputRange.hi - InputRange.lo ) .+InputRange.lo;
    outSamps = f.(samps)
    sampsMems = membership.(outFuzzys, outSamps)

    inMems = range(0,stop=1, length = Ndes)

    inLow = minimum(samps[sampsMems .!=0]); inHi = maximum(samps[sampsMems .!=0])
    inRange = interval(inLow, inHi)

    inIntervals = Interval{Float64}[]
    push!(inIntervals, inRange)

    for i = 2:Ndes

        these = findall(inMems[i-1] .<= sampsMems .<= inMems[i])

        if isempty(these); 
            push!(inIntervals, inIntervals[end]);
        else
            inLow = minimum(samps[these]); inHi = maximum(samps[these])
            inInt = interval(inLow, inHi)
            push!(inIntervals, inInt)
        end
    end

    return Fuzzy(inIntervals)
end

function invertUniv2(f :: Function, outFuzzys :: FuzzyNumber, InputRange, Nsamples = 10^4, Ndes = 200)

    samps = rand(Nsamples) .* (InputRange.hi - InputRange.lo ) .+InputRange.lo;
    outSamps = f.(samps)
    sampsMems = membership.(outFuzzys, outSamps)

    maxMem = findmax(sampsMems)
    core = interval(minimum(samps[maxMem[2]]), maximum(samps[maxMem[2]]))

    inMems = range(0,stop=1, length = Ndes)

    inLow = minimum(samps[sampsMems .!=0]); inHi = maximum(samps[sampsMems .!=0])
    inRange = interval(inLow, inHi)

    inIntervals = Interval{Float64}[]
    push!(inIntervals, inRange)

    for i = 2:Ndes

        these = findall(sampsMems .<= inMems[i])

        if isempty(these)
            push!(inIntervals, inIntervals[end]);
        else
            theseSamps = samps[these]
            leftSamps = theseSamps[theseSamps .<= core.lo]
            rightSamps = theseSamps[theseSamps .>= core.hi]

            if isempty(leftSamps); leftSamps = inIntervals[end].lo;end
            if isempty(rightSamps); rightSamps = inIntervals[end].hi;end

            push!(inIntervals, interval(maximum(leftSamps),minimum(rightSamps)))
        end
        
    end

    return Fuzzy(inIntervals)
end

function invertUnivScatter(f :: Function, outFuzzys :: FuzzyNumber, InputRange, Nsamples = 10^4, Ndes = 200)

    samps = rand(Nsamples) .* (InputRange.hi - InputRange.lo ) .+InputRange.lo;
    outSamps = f.(samps)
    sampsMems = membership.(outFuzzys, outSamps)

    return samps, sampsMems
    
end


function invertSampling(f, outData, InputRange, Nsamples = 10^4)

    Ndims = length(InputRange)

    samps = rand(Nsamples, Ndims) .* (right.(InputRange) - left.(InputRange))' .+left.(InputRange)';

    outFuzzys = [ecdf2fuzzy(outData[:,i]) for i = 1:size(outData)[2]]
    outSamps = f(samps)

    NdimsOut = length(outFuzzys)

    sampsMems = zeros(Nsamples, NdimsOut)

    for i = 1:NdimsOut
        sampsMems[:,i] = membership.(outFuzzys[i], outSamps[:,i])
    end

    joints = indepJoint(sampsMems)

    return samps, joints

end

function poss2Fuzzy1(xs, poss, Ndes = 200)

    maxMem = findmax(poss)
    core = interval(minimum(xs[maxMem[2]]), maximum(xs[maxMem[2]]))

    inMems = range(0,stop=1, length = Ndes)

    inLow = minimum(xs[poss .!=0]); inHi = maximum(xs[poss .!=0])
    inRange = interval(inLow, inHi)

    inIntervals = Interval{Float64}[]
    push!(inIntervals, inRange)

    for i = 2:Ndes

        these = findall(poss .<= inMems[i])

        if isempty(these)
            push!(inIntervals, inIntervals[end]);
        else
            theseSamps = samps[these]
            leftSamps = theseSamps[theseSamps .<= core.lo]
            rightSamps = theseSamps[theseSamps .>= core.hi]

            if isempty(leftSamps); leftSamps = inIntervals[end].lo;end
            if isempty(rightSamps); rightSamps = inIntervals[end].hi;end

            push!(inIntervals, interval(maximum(leftSamps),minimum(rightSamps)))
        end
        
    end

    return Fuzzy(inIntervals)
end


function poss2Fuzzy(xs, poss, Ndes = 200)

    maxMem = findmax(poss)
    core = interval(minimum(xs[maxMem[2]]), maximum(xs[maxMem[2]]))

    inMems = range(0,stop=1, length = Ndes)

    inLow = minimum(xs[poss .!=0]); inHi = maximum(xs[poss .!=0])
    inRange = interval(inLow, inHi)

    inIntervals = Interval{Float64}[]
    push!(inIntervals, inRange)
    
    for i = 2:Ndes
    
        these = findall(inMems[i-1] .<= poss .<= inMems[i])
    
        if isempty(these); 
            push!(inIntervals, inIntervals[end]);
        else
            inLow = minimum(xs[these]); inHi = maximum(xs[these])
            inInt = interval(inLow, inHi)
            push!(inIntervals, inInt)
        end
    end

    return Fuzzy(inIntervals)
end


inIntervals = Interval{Float64}[]
push!(inIntervals, inRange)

for i = 2:Ndes

    these = findall(inMems[i-1] .<= sampsMems .<= inMems[i])

    if isempty(these); 
        push!(inIntervals, inIntervals[end]);
    else
        inLow = minimum(samps[these]); inHi = maximum(samps[these])
        inInt = interval(inLow, inHi)
        push!(inIntervals, inInt)
    end
end

##
#   A is a Nsamples × Ndim matrix of independent possibility values
##
function indepJoint(A)

    Nsamps,Ndim = size(A)

    joint = zeros(Nsamps)

    for i =1:Nsamps
        joint[i] = 1 - (1 - minimum(A[i,:]))^Ndim
    end

    return joint
end


function f(A)
    x = A[:,1]; y =A[:,2];
    z = x.+y
    h = x.*y

    return [z h]
end