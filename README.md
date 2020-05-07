# DiffusionDefinition

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaDiffusionBayes.github.io/DiffusionDefinition.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaDiffusionBayes.github.io/DiffusionDefinition.jl/dev)
[![Build Status](https://travis-ci.com/JuliaDiffusionBayes/DiffusionDefinition.jl.svg?branch=master)](https://travis-ci.com/JuliaDiffusionBayes/DiffusionDefinition.jl)

A collection of convenient methods for defining diffusion processes. The package is designed as an integral part of a suite of packages [DiffusionBayes.jl](https://github.com/JuliaDiffusionBayes/DiffusionBayes.jl) used for Bayesian inference for diffusion processes, but it can also be used on its own to define and sample from diffusion processes.

## Minimal example
For instance, to define a Lorenz system it is enough to write
```julia
using DiffusionDefinition, StaticArrays
const DD = DiffusionDefinition

@diffusion_process Lorenz begin
    :dimensions
    process --> 3
    wiener --> 3

    :parameters
    p --> (3, Float64)
    σ --> Float64
end

DD.b(t, x, P::Lorenz) = @SVector [
    P.p1*(x[2]-x[1]),
    P.p2*x[1] - x[2] - x[1]*x[3],
    x[1]*x[2] - P.p3*x[3]
]

DD.σ(t, x, P::Lorenz) = SDiagonal(P.σ, P.σ, P.σ)
```
It is also possible to include template parameters, for instance:
```julia
@diffusion_process Lorenz{T} begin
    :dimensions
    process --> 3
    wiener --> 3

    :parameters
    p --> (3, T)
    σ --> Float64
end
```
We also provide some examples of pre-defined diffusion processes that can be
imported without having to write any code, for instance:
```
@load_diffusion :lorenz
```
and to see a list of all pre-defined examples call
```
@load_diffusion
```
We additionally provide some functionality for sampling trajectories. For instance, to sample a three-dimensional standard Brownian motion on $$[0,1]$$ use:
```julia
tt = 0.0:0.01:1.0
wiener_path = rand(Wiener(3, Float64), tt)
```
The wiener path can then be used in an Euler-Maruyama scheme to compute a trajectory under a given diffusion law:
```julia
XX = trajectory(tt, DD.ℝ{3})
P = Lorenz(28.0, 10.0, 8.0/3.0, 2.0)
DD.solve!(XX, wiener_path, P, zero(DD.ℝ{3}))
```
The two-step sampling procedure above may also be completed at once by calling:
```julia
XX = rand(P, tt, zero(DD.ℝ{3}))
```
See the [documentation](https://JuliaDiffusionBayes.github.io/DiffusionDefinition.jl/dev)
for a comprehensive overview.
