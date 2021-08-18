using PossibilisticArithmetic, Distributions, BivariateCopulas, PyPlot, Test

a = Fuzzy(1,2,3)
b = Fuzzy(1,2,4)

Cop = Pi()
Cop2 = Pi()

op = +

@time c = sigmaFuzzy(a,b,op = op, C = Cop)
@time c2 = sigmaFuzzyOld(a, b,op = +, C = Cop)

plot(c, name ="same", col = "blue")
plot(c2, name ="same", col = "red")

aDist = Uniform(1,2)
bDist = Uniform(1,2)

a2Dist = Uniform(2,3)
b2Dist = Uniform(2,3)

j1 = Cop2(aDist, bDist)
j2 = Cop2(a2Dist, b2Dist)

Nsamps = 10^5

samps = BivariateCopulas.sample(j1, Nsamps)
out = op.(samps[:,1], samps[:,2])
sort!(out)

samps2 = BivariateCopulas.sample(j2, Nsamps)
out2 = op.(samps2[:,1], samps2[:,2])
sort!(out2)

is = range(0,1,length = Nsamps)

plot(c,name = "same", col = "red")

step(out,is)
step(out2,is)

step(out,reverse(is))
step(out2,reverse(is))

probsIn1 = [sum(out .∈ cut)/Nsamps for cut in c.Membership]
probsIn2 = [sum(out2 .∈ cut)/Nsamps for cut in c.Membership]

intervalProbs = mass.(c, c.Membership)

@test all(probsIn1 .∈ intervalProbs)
@test all(probsIn2 .∈ intervalProbs)
