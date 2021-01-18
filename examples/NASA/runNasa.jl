
using DelimitedFiles


function runNasaY(a , e)
    open("aleatory.dat", "w") do io
        writedlm(io, a)
    end
    open("epistemic.dat", "w") do io
        writedlm(io, e)
    end

    try
        run(`./uqsim 1 aleatory.dat epistemic.dat`)
    finally
        return reshape(readdlm("yout"), 5001)
    end

end