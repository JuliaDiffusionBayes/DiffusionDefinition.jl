# DiffusionDefinition

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaDiffusionBayes.github.io/DiffusionDefinition.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaDiffusionBayes.github.io/DiffusionDefinition.jl/dev)
[![Build Status](https://travis-ci.com/JuliaDiffusionBayes/DiffusionDefinition.jl.svg?branch=master)](https://travis-ci.com/JuliaDiffusionBayes/DiffusionDefinition.jl)

Utility functions for defining diffusion processes in the fewest lines of code. The package is designed to work seamlessly with [BridgeSDEInference.jl](https://github.com/mmider/BridgeSDEInference.jl).

## Minimal example
For instance, to define a Lorenz system it is enough to write
```julia
@diffusion_process Lorenz begin
    :dimensions
    process --> 3
    wiener --> 3

    :parameters
    _ --> (3, Float64)
    σ --> Float64

    :additional
    constdiff --> true
end

b(t, x, P::Lorenz) = @SVector [
    P.p1*(x[2]-x[1]),
    P.p2*x[1] - x[2] - x[1]*x[3],
    x[1]*x[2] - P.p3*x[3]
]

σ(t, x, P::Lorenz) = SDiagonal(P.σ, P.σ, P.σ)
```
It is also possible to include template parameters, for instance:
```julia
@diffusion_process Lorenz{T} begin
    :dimensions
    process --> 3
    wiener --> 3

    :parameters
    _ --> (3, T)
    σ --> Float64

    :additional
    constdiff --> true
end
```
We also provide some examples of pre-defined diffusion processes that can be
imported without having to write any code with:
```
@load_diffusion :lorenz
```
To see a list of all pre-defined examples call
```
@load_diffusion
```
We additionally provide some functionality for sampling trajectories. For instance, to sample a three-dimensional standard Brownian motion use:
```julia
const DD = DiffusionDefinition
tt = collect(0.0:0.01:1.0)
wiener_path = rand(Wiener(), zero(DD.ℝ{3}), tt)
```
The wiener path can then be used in an Euler-Maruyama scheme to compute a trajectory under a given diffusion law:
```julia
XX = trajectory(tt, DD.ℝ{3})
P = Lorenz(28.0, 10.0, 8.0/3.0, 2.0)
DD.solve!(XX, wiener_path, P, zero(DD.ℝ{3}))
```
See the [documentation](https://JuliaDiffusionBayes.github.io/DiffusionDefinition.jl/stable)
for a comprehensive overview.
