# PossibilisticArithmetic.jl
![Build Status](https://github.com/Institute-for-Risk-and-Uncertainty/FuzzyArithmetic.jl/workflows/CI/badge.svg)

[Possibility Theory](https://en.wikipedia.org/wiki/Possibility_theory): fuzzy set theory applied to bounding sets of probability measures ([Imprecise Probabilities](https://en.wikipedia.org/wiki/Imprecise_probability)). This package gives a rigorous arithmetic between fuzzy numbers. Outer approximations and arithmetic between fuzzy sets are based on [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl).

### In development. Package can currently:

* Get interval bounds on probabilities
* Construct robust outer approximations of membership functions with different step sizes
* Construct fuzzy numbers with a real or interval core
* Perform dependent arithmeric (using [copulas](https://github.com/AnderGray/BivariateCopulas.jl)) between fuzzy numbers
    * Levelwise arithmetic (standard fuzzy arithmetic)
    * Arithmetic with a precise copula defined (sigma convolutions)
    * Arithmetic with a copula's bound defined (tau convolutions)
* Unary operators: `-, sin, cos, tan, exp, log`
* Interacts with [ProbabilityBoundsAnalysis.jl](https://github.com/AnderGray/ProbabilityBoundsAnalysis.jl) and [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl)
* Conversions between fuzzy numbers, [p-boxes](https://en.wikipedia.org/wiki/Probability_box), and more general [random sets](https://en.wikipedia.org/wiki/Dempster–Shafer_theory)
* Plot fuzzy numbers (requires PyPlot)

### Soon:
* Multimodal possibility distributions (non-convex fuzzy sets) using [IntervalUnionArithmetic.jl](https://github.com/AnderGray/IntervalUnionArithmetic.jl)
* Parameteric distribution constructors
* Documentation and testing

Resources
---

* Fuzzy & Possibility Interest Group (FnPIG) [google site](https://sites.google.com/site/fuzzypossrisk/) (access granted on request)
* Society for Imprecise Probabilities: Theories and Applications [website](https://www.sipta.org)
* [ISIPTA21 paper](https://proceedings.mlr.press/v147/gray21a.html)
* [ISIPTA21 poster](https://www.researchgate.net/publication/353220811_Poster_Dependent_Possibilistic_Arithmetic_using_Copulas)


Installation
---

`PossibilisticArithmtic.jl` is a registered Julia package, and can be installed using the Julia package manager:
```julia
julia> ]
(v1.0) pkg> add PossibilisticArithmetic
julia> using PossibilisticArithmetic
```

or for the latest version:

```julia
julia> ]
(v1.0) pkg> add https://github.com/AnderGray/PossibilisticArithmetic.jl#master
julia> using PossibilisticArithmetic
```

Use
---

```julia
julia> a = Fuzzy(0, 1, 2)
Fuzzy: 	 ~ ( Range=[0, 2], Core=1.0 )

julia> a = Fuzzy(0, 1, 2, steps = 200) # Give number of steps, default = 200
Fuzzy: 	 ~ ( Range=[0, 2], Core=1.0 )

julia> cut(a,0.2)       # get an Alpha-cut, returns an interval
[0.195979, 1.80403]

julia> b = Fuzzy(0, 1, 1.5, 2)
Fuzzy: 	 ~ ( Range=[0, 2], Core=[1, 1.5] )

julia> c = a * b
Fuzzy: 	 ~ ( Range=[0, 4], Core=[1, 1.5] )

julia> mass(c, interval(3, 4) )   # Bound probability measure in some set
[0, 0.351759]

julia> mass(c, 3, 4)   # Same as above
[0, 0.351759]

julia> c(3, 4)         # Can also query object directly
[0, 0.351759]

julia> plot(a)
julia> plot(b)
julia> plot(c)
```
![Fuzzy with real core](https://i.imgur.com/7ZYbTyR.png)
![Fuzzy with interval core](https://i.imgur.com/h8h3u7c.png)
![Levelwise fuzzy product](https://i.imgur.com/pq4djBT.png)
