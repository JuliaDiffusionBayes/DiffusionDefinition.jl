# Loading diffusions
We provide definitions of many standard diffusion processes. To see a list of all available processes call
```julia
@load_diffusion
```
See the [Examples](/examples/lotka_volterra) for more detailed descriptions. To load a particular diffusion simply write: `@load_diffusion DiffusionName`. For instance:
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
