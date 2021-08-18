using FuzzyArithmetic, BivariateCopulas, ProbabilityBoundsAnalysis, LaTeXStrings


a = Fuzzy(1,2,3)
b = Fuzzy(1,2,3)

plot(a, name ="afuz", fontsize= 24)
PyPlot.title(L"$\pi_X, \pi_Y$", fontsize = 60)
PyPlot.xlabel(L"$X,Y$", fontsize = 32)
PyPlot.ylabel(L"$\alpha$", fontsize = 32)
savefig("newFig4/aFuz.pdf")

aBox = makepbox(a)

ProbabilityBoundsAnalysis.plot(aBox, name = "abox", fontsize = 24)
PyPlot.title(L"$\mathbf{F}_{X}, \mathbf{F}_Y$", fontsize = 60)
PyPlot.xlabel(L"$X,Y$", fontsize = 32)
PyPlot.ylabel(L"$CDF$", fontsize = 32)
savefig("newFig4/abox.pdf")

op = +

int = op.(a.Membership[1],b.Membership[end])

c1 = Fuzzy([int])


cors = [-1, -0.8, -0.3 , 0, 0.3, 0.8, 1]
cols = ["red", "yellow", "green", "blue", "orange", "black", "purple"]

cops = Gaussian.(cors)

# fuzzy sigma

plot(c1, name = "sigmaFuzzy", col = "red", fontsize = 24)


for i =2:length(cors)
    C1 = cops[i]
    this = sigmaFuzzy(a,b, op = op, C=C1)
    plot(this, name ="sigmaFuzzy", col = cols[i],fontsize = 24)
end

PyPlot.title(L"$\pi_Z$", fontsize = 60)
PyPlot.xlabel(L"$Z$", fontsize = 32)
PyPlot.ylabel(L"$\alpha$", fontsize = 32)
savefig("newFig4/sigmaFuzzy.pdf")

thisbox = ProbabilityBoundsAnalysis.makepbox([int])
ProbabilityBoundsAnalysis.plot(thisbox,  name = "sigmaFuzzyBox", col = "red", fontsize = 24)

for i =2:length(cors)
    C1 = cops[i]
    this = sigmaFuzzy(a,b, op = op, C=C1)

    thisbox = makepbox(this)
    ProbabilityBoundsAnalysis.plot(thisbox,name= "sigmaFuzzyBox", col = cols[i],fontsize = 24)

end
PyPlot.title(L"$\mathbf{F}_Z$", fontsize = 60)
PyPlot.xlabel(L"$Z$", fontsize = 32)
PyPlot.ylabel(L"$CDF$", fontsize = 32)
savefig("newFig4/sigmaBox.pdf")

# Tau sigma


for i =1:length(cors)
    C1 = cops[i]
    this = tauFuzzy(a,b, op = op, C=C1)
    plot(this, name ="tauFuzzy", col = cols[i],fontsize = 24)

end
PyPlot.title(L"$\pi_Z$", fontsize = 60)
PyPlot.xlabel(L"$Z$", fontsize = 32)
PyPlot.ylabel(L"$\alpha$", fontsize = 32)
savefig("newFig4/tauFuzz.pdf")

for i =1:length(cors)
    C1 = cops[i]
    this = tauFuzzy(a,b, op = op, C=C1)

    thisbox = makepbox(this)
    ProbabilityBoundsAnalysis.plot(thisbox,name= "tauFuzzyBox", col = cols[i],fontsize = 24)

end
PyPlot.title(L"$\mathbf{F}_Z$", fontsize = 60)
PyPlot.xlabel(L"$Z$", fontsize = 32)
PyPlot.ylabel(L"$CDF$", fontsize = 32)
savefig("newFig4/tauBox.pdf")

# pba calculations

abox = makepbox(a);
bbox = makepbox(b);

ProbabilityBoundsAnalysis.plot(abox, name = "abox", fontsize = 24)
ProbabilityBoundsAnalysis.plot(bbox, name = "bbox", fontsize = 24)

for i =1:length(cors)
    C1 = ProbabilityBoundsAnalysis.GauCopula(cors[i])
    this = ProbabilityBoundsAnalysis.tauRho(abox,bbox, op = op, C=C1)
    ProbabilityBoundsAnalysis.plot(this, name ="tauPbox", col = cols[i],fontsize = 24)

    that = sigma(abox,bbox, op = op, C=C1)

    ProbabilityBoundsAnalysis.plot(that, name="sigmaBox", col = cols[i],fontsize = 24)

end
