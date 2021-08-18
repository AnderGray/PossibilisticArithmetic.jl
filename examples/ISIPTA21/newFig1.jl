using PossibilisticArithmetic, BivariateCopulas, LaTeXStrings

a = Fuzzy(1,2,3)
b = Fuzzy(1,2,3)

cors = [-1, -0.8, -0.3 , 0, 0.3, 0.8, 1]
cols = ["red", "yellow", "green", "blue", "orange", "black", "purple"]

cops = Gaussian.(cors)

op = *

for i =1:length(cors)
    this = sigmaFuzzy(a,b, op = op, C=cops[i])
    plot(this, name ="X * Y", col = cols[i],fontsize = 24)
end
PyPlot.title(L"$\pi_X * \pi_Y$", fontsize = 32)
PyPlot.xlabel(L"$Z$", fontsize = 24)
PyPlot.ylabel(L"$\alpha$", fontsize = 24)
savefig("sigmaNew/X * Y.pdf")


cols = reverse(cols)
op = -
for i =2:length(cors)
    this = sigmaFuzzy(a,b, op = op, C=cops[i])
    plot(this, name ="X - Y ", col = cols[i-1],fontsize = 24)
end

int = op(a.Membership[1], b.Membership[end])

this = Fuzzy([int])
plot(this, name ="X - Y ", col = "purple",fontsize = 24)

PyPlot.title(L"$\pi_X - \pi_Y$", fontsize = 32)
PyPlot.xlabel(L"$Z$", fontsize = 24)
PyPlot.ylabel(L"$\alpha$", fontsize = 24)
savefig("sigmaNew/X - Y.pdf")


op = /

for i =1:length(cors)
    this = sigmaFuzzy(a,b, op = op, C=cops[i])
    plot(this, name ="X / Y", col = cols[i],fontsize = 24)
end

PyPlot.title(L"$\pi_X \div \pi_Y$", fontsize = 32)
PyPlot.xlabel(L"$Z$", fontsize = 24)
PyPlot.ylabel(L"$\alpha$", fontsize = 24)
savefig("sigmaNew/X div Y.pdf")

#=
op = min

for i =1:length(cors)
    this = sigmaFuzzy(a,b, op = op, C=cops[i])
    plot(this, name ="min(X, Y)", col = cols[i],fontsize = 24)
end

PyPlot.title("min(X, Y)", fontsize = 24)
savefig("sigmaNew/min(X, Y)")


op = max

for i =1:length(cors)
    this = sigmaFuzzy(a,b, op = op, C=cops[i])
    plot(this, name ="max(X, Y)", col = cols[i],fontsize = 24)
end

PyPlot.title("max(X, Y)", fontsize = 24)
savefig("sigmaNew/max(X, Y)")


thisop(x,y) = sin(x) + y^2

for i =1:length(cors)
    this = sigmaFuzzy(a,b, op = thisop, C=cops[i])
    plot(this, name ="sin(x) + y^2", col = cols[i],fontsize = 24)
end

PyPlot.title("sin(x) + y^2", fontsize = 24)
savefig("sigmaNew/sin(x) + y^2")


thisop(x,y) = sin(x) + cos(y)

for i =1:length(cors)
    this = sigmaFuzzy(a,b, op = thisop, C=cops[i])
    plot(this, name ="sin(x) + cos(y)", col = cols[i],fontsize = 24)
end

PyPlot.title("sin(x) + cos(y)", fontsize = 24)
savefig("sigmaNew/sin(x) + cos(y)")


thisop(x,y) = exp(x) * y^2

for i =1:length(cors)
    this = sigmaFuzzy(a,b, op = thisop, C=cops[i])
    plot(this, name ="exp(x) * y^2", col = cols[i],fontsize = 24)
end

PyPlot.title("exp(x) * y^2", fontsize = 24)
savefig("sigmaNew/exp(x) * y^2")

=#
