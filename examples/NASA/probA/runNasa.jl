

using Distributed
@everywhere using DelimitedFiles
@everywhere using Base

function runNasaY(a, e)

    Nsamples = size(a,2)

    open("aleatory.dat", "w") do io
        writedlm(io, a)
    end
    open("epistemic.dat", "w") do io
        writedlm(io, e)
    end

    run(`./uqsim 1 aleatory.dat epistemic.dat`)

    return reshape(readdlm("yout"), Nsamples,5001)'
end

function runNasaY(Nsamples)

    a = rand(5,Nsamples) .*2
    e = rand(4,Nsamples) .*2

    Nsamples = size(a,2)

    open("aleatory.dat", "w") do io
        writedlm(io, a)
    end
    open("epistemic.dat", "w") do io
        writedlm(io, e)
    end

    run(`./uqsim 1 aleatory.dat epistemic.dat`)

    return reshape(readdlm("yout"), Nsamples, 5001)'
end

function runNasaYFFT(Nsamples)

    a = rand(5,Nsamples) .*2
    e = rand(4,Nsamples) .*2

    Nsamples = size(a,2)

    open("aleatory.dat", "w") do io
        writedlm(io, a)
    end
    open("epistemic.dat", "w") do io
        writedlm(io, e)
    end

    run(`./uqsim 1 aleatory.dat epistemic.dat`)
    ys = reshape(readdlm("yout"),Nsamples, 5001)'

    return abs.(fft(ys))[1:50,:]
end


function runNasaYFFT2(A)

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
    ys = reshape(readdlm("yout"),5001,Nsamples)'

    return abs.(fft(ys))[1:50,:]
end


@everywhere function runNasaYDir(a, e)

    dirname = "Nasa_1"
    i=1
    while isdir(dirname)
        i = i+1
        dirname = "Nasa_$i"
    end
    mkpath(dirname)
    cd(dirname)
    run(`cp ../uqsim .`)
    
    Nsamples = size(a,2)

    open("aleatory.dat", "w") do io
        writedlm(io, a)
    end
    open("epistemic.dat", "w") do io
        writedlm(io, e)
    end

    run(`./uqsim 1 aleatory.dat epistemic.dat`)
    cd("..")
    return reshape(readdlm("$dirname/yout"), Nsamples, 5001)';
end

@everywhere function runNasaYDir2(Nsamples, id)

    dirname = "Nasa_$id"
    cd(dirname)
    
    
    a = rand(5,Nsamples) .*2
    e = rand(4,Nsamples) .*2

    open("aleatory.dat", "w") do io
        writedlm(io, a)
    end
    open("epistemic.dat", "w") do io
        writedlm(io, e)
    end

    run(`./uqsim 1 aleatory.dat epistemic.dat`)
    cd("..")
    return reshape(readdlm("$dirname/yout"), 5001, Nsamples)
end



function runNasaYPar(Nsamples, Ncores)
    
    [isdir("Nasa_$i") ? println("Nasa_$i exists") : mkdir("Nasa_$i")  for i =1:Ncores]
    [run(`cp uqsim Nasa_$i`) for i =1:Ncores]

    N = Integer(Nsamples/Ncores)

    jobs = [Distributed.@spawnat :any runNasaYDir2(N,i) for i = 1:Ncores]

    #outs = [fetch(j) for j in jobs]

    return [fetch(j) for j in jobs]
end