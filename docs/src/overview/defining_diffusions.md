# Defining diffusion
The main utility macro introduced in this package is `@diffusion_process`. It lets us define a struct characterising a diffusion process in a minimal amount of written code.

There are two (optionally three, third is not used now) parts expected by the macro: the name of the diffusion and the recipe for defining a struct (the third one, not used we denote as POSTVARIABLES):
```julia
@diffusion_process NAME begin
  RECIPE
end POSTVARIABLES
```
snippet of code above creates a struct named `NAME` according to specifications listed in `RECIPE`.

## Customisation of a struct
The `RECIPE` can contain information of four different types:
- Specification of `:dimensions`
- Specification of `:parameters` (their names and datatypes)
- Specification of `:conjugate` updates
- `:additional` information
Each type needs to be announced to `julia` by starting the list with the corresponding `Symbol` (or `QuoteNode`).
### Example
It's best to look at a simple example. Consider definition of a Lorenz system:
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
```
