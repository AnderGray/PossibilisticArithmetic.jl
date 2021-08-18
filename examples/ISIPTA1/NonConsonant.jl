using FuzzyArithmetic, BivariateCopulas,LaTeXStrings

function sigmaFuzzyNC(x::Fuzzy, y::Fuzzy; op = +, C = Pi())

    #if op = -; return sigmaFuzzy(x, -y, op=+, C = rotate(C)); end
    #if op = /; return sigmaFuzzy(x, 1/y, op=/, C = rotate(C)); end

    if C == W; return levelwiseOpp(x, y, op = op);end
    if C == M; return levelwise(x, y, op = op);end

    zNum = max(length(x.Membership), length(y.Membership))
    masses, cartProd = mobiusTransform2D(x, y, C)       # Get carteesian prod and masses from Möbius

    zs = [op(ints[1], ints[2]) for ints in cartProd]     # Evaluate carteesian product with interval arithm

    zs = sort(zs,lt = ⊆)

    return Fuzzy( reverse(zs))
end

function sigmaFuzzySlow(x::Fuzzy, y::Fuzzy; op = +, C = Pi())

    #if op = -; return sigmaFuzzy(x, -y, op=+, C = rotate(C)); end
    #if op = /; return sigmaFuzzy(x, 1/y, op=/, C = rotate(C)); end

    if C == W; return levelwiseOpp(x, y, op = op);end
    if C == M; return levelwise(x, y, op = op);end

    zNum = max(length(x.Membership), length(y.Membership))
    masses, cartProd = mobiusTransform2D(x, y, C)       # Get carteesian prod and masses from Möbius

    zs = [op(ints[1], ints[2]) for ints in cartProd]     # Evaluate carteesian product with interval arithm

    return FuzzyArithmetic.DSS2FuzzySlow(zs, masses,steps= zNum)
end


a = Fuzzy(1, 2, 3, steps = 100)
b = Fuzzy(1, 2, 10, steps = 100)

op = *

c = sigmaFuzzyNC(a,b, op = op)
c2 = sigmaFuzzySlow(a,b, op = op)
c3 = sigmaFuzzy(a,b,op=op)

cuts = c2.Membership

plot(c, fontsize = 24)
PyPlot.xlabel("Z", fontsize = 24)
title(L"$\mathrm{Bel}_Z$", fontsize = 32)


plot(c2, fontsize = 24)
PyPlot.xlabel("Z", fontsize = 24)
title(L"$\pi_Z$", fontsize = 32)


m1 = mass.(c,cuts)
m2 = mass.(c2,cuts)

m3 = mass.(c, c.Membership)
m4 = mass.(c2, c.Membership)

m1 .⊆ m2

all(m1 .⊆ m2)
all(m3 .⊆ m4)


findall( .~(m1 .⊆ m2))
