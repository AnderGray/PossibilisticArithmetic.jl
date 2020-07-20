# FuzzyArithmetic.jl

A julia package for performing fuzzy arithmetic.

### In development. Package can currently:

* Construct fuzzy numbers with a real or interval core
* Have different discretization of the membership function
* Perform alpha cuts
* Perform levelwise arithmeric between fuzzy numbers, intervals and scalars
* Get bounds on probability mass in an interval
* Plot fuzzy numbers

## Install:

```BASH
# Download from git
git clone https://github.com/Institute-for-Risk-and-Uncertainty/FuzzyArithmetic.jl.git  

# Enter package
cd FuzzyArithmetic.jl       

# Begin julia
julia                       
```
```julia
julia> include("src/FuzzyNumbers.jl")
```

## Use:

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
