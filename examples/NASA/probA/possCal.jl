
@everywhere include("readData.jl")
@everywhere include("runNasa.jl")
@everywhere include("../../../src/FuzzyNumbers.jl")


@everywhere function runNasaYFFT2(A, id)

    dirname = "Nasa_$id"
    cd(dirname)

    a = A[1:5,:]
    e = A[6:9,:]

    Nsamples = size(a,2)

    open("aleatory.dat", "w") do io
        writedlm(io, a)
    end
    open("epistemic.dat", "w") do io
        writedlm(io, e)
    end

    run(`./uqsim 1 aleatory.dat epistemic.dat`)
    ys = reshape(readdlm("yout"), 5001,Nsamples)
    return abs.(fft(ys))[1:50,:]
end

@everywhere function invertSampling(f, outData, InputRange, id, Nsamples = 10^4)

    Ndims = length(InputRange)

    samps = rand(Nsamples, Ndims) .* (right.(InputRange) - left.(InputRange))' .+left.(InputRange)';

    outFuzzys = [ecdf2fuzzy(outData[i,:]) for i = 1:size(outData)[1]]
    outSamps = f(samps', id)

    NdimsOut = length(outFuzzys)

    sampsMems = zeros(Nsamples, NdimsOut)

    for i = 1:NdimsOut
        sampsMems[:,i] = membership.(outFuzzys[i], outSamps[i,:])
    end

    joints = indepJoint(sampsMems)

    open("possibilities.dat", "w") do io
        writedlm(io, joints)
    end

end
@everywhere function indepJoint(A)

    Nsamps,Ndim = size(A)

    joint = zeros(Nsamps)

    for i =1:Nsamps
        joint[i] = 1 - (1 - minimum(A[i,:]))^Ndim
    end

    return joint
end


yData = getMatData()

yFFT = abs.(fft(yData))[1:50,:]

inRange = interval(0,2)×interval(0,2)×interval(0,2)×interval(0,2)×interval(0,2)×interval(0,2)×interval(0,2)×interval(0,2)×interval(0,2)

Ncores = 10

jobs = [Distributed.@spawnat :any invertSampling(runNasaYFFT2, yFFT, inRange, i, 10^4) for i in 1:Ncores]

s = [fetch(j) for j in jobs]