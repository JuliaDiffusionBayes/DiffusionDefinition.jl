# DiffusionDefinition

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mmider.github.io/DiffusionDefinition.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mmider.github.io/DiffusionDefinition.jl/dev)
[![Build Status](https://travis-ci.com/mmider/DiffusionDefinition.jl.svg?branch=master)](https://travis-ci.com/mmider/DiffusionDefinition.jl)

[PACKAGE UNDER DEVELOPMENT]

A set of utility functions for defining diffusion processes with minimal code. Designed to work seamlessly with [BridgeSDEInference.jl](https://github.com/mmider/BridgeSDEInference.jl).

## Minimal example
For instance, to define a Lorenz system it is enough to write
```julia
@diffusion_process Lorenz begin
    :dimensions
    #---------
    process --> 3
    wiener --> 3

    :parameters
    #---------
    _ --> (3, Float64)
    σ --> Float64

    :additional
    #---------
    constdiff --> true
end

function b(t, x, P::Lorenz)
    @SVector [
        P.p1*(x[2]-x[1]),
        P.p2*x[1] - x[2] - x[1]*x[3],
        x[1]*x[2] - P.p3*x[3]
    ]
end

function σ(t, x, P::Lorenz)
    @SMatrix [
        P.σ 0.0 0.0;
        0.0 P.σ 0.0;
        0.0 0.0 P.σ
    ]
end
```
See the [documentation](https://mmider.github.io/DiffusionDefinition.jl/stable) for a comprehensive overview.