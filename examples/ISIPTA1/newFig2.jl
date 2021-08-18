using FuzzyArithmetic, BivariateCopulas, LaTeXStrings


a = Fuzzy(1,2,3)
b = Fuzzy(1,2,3)

cors = [-1, -0.8, -0.3 , 0, 0.3, 0.8, 1]
cols = ["red", "yellow", "green", "blue", "orange", "black", "purple"]

cops = Gaussian.(cors)


op = *

cors = [-1, -0.8, -0.3 , 0, 0.3, 0.8, 1]
cols = ["red", "yellow", "green", "blue", "orange", "black", "purple"]

for i =1:length(cors)
    this = tauFuzzy(a,b, op = op, C=cops[i])
    plot(this, name ="X * Y_τ", col = cols[i],fontsize = 24)
end
PyPlot.xlabel(L"$Z$", fontsize= 24)
PyPlot.ylabel(L"$\alpha$", fontsize= 24)
PyPlot.title(L"\pi_X * \pi_Y", fontsize = 32)
savefig("tauNew/X * Y_τ.pdf")


cols = reverse(cols);
op = -

for i =1:length(cors)
    this = tauFuzzy(a,b, op = op, C=cops[i])
    plot(this, name ="X - Y_τ", col = cols[i],fontsize = 24)
end
PyPlot.xlabel(L"$Z$", fontsize= 24)
PyPlot.ylabel(L"$\alpha$", fontsize= 24)
PyPlot.title(L"\pi_X - \pi_Y", fontsize = 32)
savefig("tauNew/X - Y_τ.pdf")

op = /

for i =1:length(cors)
    this = tauFuzzy(a,b, op = op, C=cops[i])
    plot(this, name ="X / Y_τ", col = cols[i],fontsize = 24)
end

PyPlot.xlabel(L"$Z$", fontsize= 24)
PyPlot.ylabel(L"$\alpha$", fontsize= 24)
PyPlot.title(L"\pi_X \div \pi_Y", fontsize = 32)
savefig("tauNew/X div Y_τ.pdf")

op = min

for i =1:length(cors)
    C1 = Gaussian(cors[i])
    this = tauFuzzy(a,b, op = op, C=cops[i])
    plot(this, name ="min(X, Y)_τ", col = cols[i],fontsize = 24)
end

PyPlot.title("min(X, Y)", fontsize = 24)
savefig("tauNew/min(X, Y)_τ")


op = max

for i =1:length(cors)
    this = tauFuzzy(a,b, op = op, C=cops[i])
    plot(this, name ="max(X, Y)_τ", col = cols[i],fontsize = 24)
end

PyPlot.title("max(X, Y)", fontsize = 24)
savefig("tauNew/max(X, Y)_τ")



thisop(x,y) = sin(x) + y^2

for i =1:length(cors)
    this = tauFuzzy(a,b, op = thisop, C=cops[i])
    plot(this, name ="sin(x) + y^2_τ", col = cols[i],fontsize = 24)
end

PyPlot.title("sin(x) + y^2", fontsize = 24)
savefig("tauNew/sin(x) + y^2_τ")




thisop(x,y) = sin(x) + cos(y)

for i =1:length(cors)
    this = tauFuzzy(a,b, op = thisop, C=cops[i])
    plot(this, name ="sin(x) + cos(y)_τ", col = cols[i],fontsize = 24)
end

PyPlot.title("sin(x) + cos(y)", fontsize = 24)
savefig("tauNew/sin(x) + cos(y)_τ")



thisop(x,y) = exp(x) * y^2

for i =1:length(cors)
    this = tauFuzzy(a,b, op = thisop, C=cops[i])
    plot(this, name ="exp(x) * y^2_τ", col = cols[i],fontsize = 24)
end

PyPlot.title("exp(x) * y^2", fontsize = 24)
savefig("tauNew/exp(x) * y^2_τ")
