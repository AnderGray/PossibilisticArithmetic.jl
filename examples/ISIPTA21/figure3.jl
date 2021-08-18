using PossibilisticArithmetic, LaTeXStrings, PyPlot, BivariateCopulas


a = Fuzzy(1,2,3)
b = Fuzzy(1,2,3)

op = +
#cors = [-1, -0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.8, 1]
#cols = ["red", "blue", "green", "yellow", "orange", "black", "purple", "red", "blue"]

#cors = [-0.8, -0.3 , 0, 0.3, 0.8, 1]
#cols = ["yellow", "green", "blue", "orange", "black", "purple"]

cors = [-1, -0.8, -0.3 , 0, 0.3, 0.8, 1]
cols = ["red", "yellow", "green", "blue", "orange", "black", "purple"]

#plot(c1,name = "same", col = "red", fontsize = 24)

for i =1:length(cors)
    C1 = Gaussian(cors[i])
    this = tauFuzzy(a,b, op = op, C=C1)
    #f = supT(a, b, op = op, T = C1)
    plot(this, name ="same", col = cols[i],fontsize = 24)
end
PyPlot.xlabel(L"$Z$", fontsize= 24)
PyPlot.ylabel(L"$\alpha$", fontsize= 24)
PyPlot.title(L"\pi_X + \pi_Y", fontsize = 32)
savefig("figure3.pdf")
