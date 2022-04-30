# NCDist.jl

## About

NCDist.jl is a Julia package to compute distributions as well as Brown measures of (self-adjoint) rational expressions and operator-valued elements given the distributions of the variables.
It provides implementations of the [subordination algorithm](https://arxiv.org/abs/1303.3196) and a [specialized algorithm](https://arxiv.org/abs/math/0703510) for semi-circle distributions.

## Installation

In the package manager (accessible in the REPL via `]`):

```julia
(@v1.7) pkg> add git@github.com:johannes-hoffmann/NCDist.jl.git
```

Now leave the package manager (via backspace) and load the package:

```julia
julia> using NCDist
```

## Usage

After loading the package, you can find detailed usage instructions including information on the parameters for all of the functions via the help functionality of the REPL by typing a question mark `?` followed by the name of the function.

### Provided functions

To compute distributions for arbitrary variable distributions:

* `distribution_rational`
* `distribution_operator`

To compute distributions for semi-circle variable distributions:

* `distribution_rational_semicircle`
* `distribution_operator_semicircle`

To compute Brown measures for arbitrary variable distributions at a collection of points:

* `brown_rational`
* `brown_operator`

To compute Brown measures for arbitrary variable distributions at points on a grid:

* `brown_rational_grid`
* `brown_operator_grid`

To compute Brown measures for semi-circle variable distributions at a collection of points:

* `brown_rational_semicircle`
* `brown_operator_semicircle`

To compute Brown measures for semi-circle variable distributions at points on a grid:

* `brown_rational_grid_semicircle`
* `brown_operator_grid_semicircle`

### Inputs

Most of the functions in this package use inputs of the following types:

#### Rational expressions

You need to provide the rational expression or operator-valued element in the form of a linearization `L`, which in this package is an `Array{Complex{Float64}, 3}`, where `L[:, :, k]` is the matrix corresponding to the variable `k - 1` (`L[:, :, 1]` corresponds to the constant part).

#### Variable distributions

Furthermore, the distributions of the variables `var_dists` are given as a `Vector{<:Function}`, whose entries are their Cauchy transforms (on the upper half-plane) . For convenience, several common cases are provided.

* `semicircle`
* `poisson`
* `bernoulli`
* `symmetric_bernoulli`
* `cauchy`
* `arcsine`

## Example

Create a linearization:

```julia
julia> L = zeros(Complex{Float64}, 2, 2, 3);
julia> L[2, 2, 2] = 1;
julia> L[1, 2, 3] = 1;
julia> L[2, 1, 3] = 1;
```

Define the variable distributions (the first variable has a semi-circle distribution with mean `a=1` and variance `t=1`, the second variable a free Poisson distribution with default parameters):

```julia
julia> var_dists = [semicircle(a=0, t=1), poisson()]
```

Compute the distribution at points `-0.01`, `0.0`, and `0.01`:

```julia
julia> distribution_rational(L, var_dists, [-0.01, 0.0, 0.01])
3-element Vector{Float64}:
 -2.41726950654163e-6
 -0.6267424782633638
 -2.41726947793608e-6
```