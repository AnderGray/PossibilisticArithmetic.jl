using PossibilisticArithmetic

a = FuzzyNumber(1,2,3)
b = FuzzyNumber(1,2,3)

abox = makepbox(a)
bbox = makepbox(b)

c = tauFuzzy(a,b)

cFbox = makepbox(c)

cbox = abox+bbox

cboxTau = tauRho(abox,bbox, op = +, C = πCop())

cSig = sigmaFuzzy(a,b)

cSigBox = makepbox(cSig)

ProbabilityBoundsAnalysis.plot(cFbox, name = "same", col = "red", fontsize = 24)
ProbabilityBoundsAnalysis.plot(cbox, name = "same", col = "blue",fontsize = 24)
#ProbabilityBoundsAnalysis.plot(cboxTau, name = "same", col = "black")
#ProbabilityBoundsAnalysis.plot(cSigBox, name = "same", col = "green")


cboxSigFuz = makefuzzy(cbox)

plot(cboxSigFuz,name = "same2", col = "red")
plot(c,name = "same2", col = "blue")
