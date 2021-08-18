# PossibilisticArithmetic.jl
![Build Status](https://github.com/Institute-for-Risk-and-Uncertainty/FuzzyArithmetic.jl/workflows/CI/badge.svg)

[Possibility Theory](https://en.wikipedia.org/wiki/Possibility_theory): fuzzy set theory applied to bounding sets of probability measures ([Imprecise Probabilities](https://en.wikipedia.org/wiki/Imprecise_probability)). This package gives a rigorous arithetic between possibility distributions.

### In development. Package can currently:

* Construct fuzzy numbers with a real or interval core
* Have different discretization of the membership function
* Get interval bounds on probabilities
* Perform dependent arithmeric (using [copulas](https://github.com/AnderGray/BivariateCopulas.jl)) between fuzzy numbers, intervals and scalars
* Unary operators: `-, sin, cos, tan, exp, log`
* Interacts with [ProbabilityBoundsAnalysis.jl](https://github.com/AnderGray/ProbabilityBoundsAnalysis.jl)
* Plot fuzzy numbers


Installation
---

`PossibilisticArithmtic.jl` is a not yet a registered Julia package, and so the latest release can be installed using the Julia package manager:

```julia
julia> ]
(v1.0) pkg> add https://github.com/AnderGray/FuzzyArithmetic.jl
```

Use
---

```julia
julia> a = FuzzyNumber(0, 1, 2)
Fuzzy: 	 ~ ( Range=[0, 2], Core=1.0 )

julia> cut(a,0.2)
[0.195979, 1.80403]

julia> b = FuzzyNumber(0, 1, 1.5, 2)
Fuzzy: 	 ~ ( Range=[0, 2], Core=[1, 1.5] )

julia> c = a * b
Fuzzy: 	 ~ ( Range=[0, 4], Core=[1, 1.5] )

julia> mass(c, 3, 4)   # Get probability mass between two values
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
