# [Installation](@id get_started)
The package is not registered yet. To install it write:
```julia
] add https://github.com/JuliaDiffusionBayes/DiffusionDefinition.jl
```
# Usage
!!! warning "TODO"
    change for something simpler

## Defining the process
To define a diffusion law use a macro `@diffusion_definition`:
```julia
using DiffusionDefinition
const DD = DiffusionDefinition

@diffusion_process Lorenz begin
    :dimensions
    process --> 3
    wiener --> 3

    :parameters
    p --> (3, Float64)
    σ --> Float64

    :additional
    constdiff --> true
    sparsediff --> true
end
```
This will define a struct `Lorenz` and announce to Julia that it represents a three-dimensional diffusion process driven by a three-dimensional Brownian motion and that it depends on 4 parameters of type `Float64` each. It will also say something useful about the volatility coefficient.

To complete characterization of a diffusion law we define the drift and diffusion coefficients:
```julia
function DD.b(t, x, P::Lorenz)
    @SVector [
        P.p1*(x[2]-x[1]),
        P.p2*x[1] - x[2] - x[1]*x[3],
        x[1]*x[2] - P.p3*x[3]
    ]
end

DD.σ(t, x, P::Lorenz) = SDiagonal(P.σ, P.σ, P.σ)
```
## Sampling the process
To sample the process we use a function `rand`
```julia
tt = 0.0:0.001:1.0
y1 = rand(DD.ℝ{3})
P = Lorenz(10.0, 8.0/3.0, 28.0, 0.2)
X = rand(P, tt, y1)
```
## Plotting the process
You may use your favourite plotting for backend and use function `plot` on the output of `rand`:
```julia
using Plots
gr()
plot(X, Val(:vs_time))
```
