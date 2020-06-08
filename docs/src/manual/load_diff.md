# Loading diffusions
******
Many standard diffusion processes have already been defined by us in this package via `@diffusion_process` and you may simply load them and start using them immediately without having to define them yourself. To see a list of all available processes call
```julia
@load_diffusion
```
See the [Predefined processes](@ref example_lotka_volterra) for more detailed descriptions. To load a particular diffusion simply write: `@load_diffusion DiffusionName`. For instance:
```julia
@load_diffusion LotkaVolterraAux
```
Then the process can be instantiated:
```julia
α, β, γ, δ, σ1, σ2 = 2.0/3.0, 4.0/3.0, 1.0, 1.0, 0.2, 0.3
lotka_volterra = LotkaVolterraAux(
    α, β, γ, δ, σ1, σ2,
    0.0, 1.0, zero(DD.ℝ{2}), zero(DD.ℝ{2})
)
```
A very small subgroup of predefined diffusions have a variable dimension parameter. For these we provide a `@load_variable_diffusion` macro, that in addition to a diffusion name accepts also a dimension of the process and dimension of the wiener process. For instance
```julia
@load_variable_diffusion Lorenz96 2^10
# if it made sense an additional 2nd number would specify the dimension of a Wiener process
```
