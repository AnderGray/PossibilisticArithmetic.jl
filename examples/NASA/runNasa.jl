
using DelimitedFiles


function runNasaY(a , e)
    open("aleatory.dat", "w") do io
        writedlm(io, a)
    end
    open("epistemic.dat", "w") do io
        writedlm(io, e)
    end

    run(`./uqsim 1 aleatory.dat epistemic.dat design.dat`)
    
    return reshape(readdlm("yout"), 5001)

end