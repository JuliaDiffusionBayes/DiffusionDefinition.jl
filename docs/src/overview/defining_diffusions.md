# Defining diffusion
The main utility macro introduced in this package is `@diffusion_process`. It facilitates very concise definitions of structs characterizing diffusion processes.

There are two (optionally three, third is not used now) parts expected by the macro:
- the name of the diffusion (which may contain template parameters in the curly brackets)
- and the recipe for defining a struct
- (the third one, not used for now, we denote as POSTVARIABLES):
```julia
@diffusion_process NAME{TEMPLATE_PARAMETERS} begin
  RECIPE
end POSTVARIABLES
```
the snippet of code above creates a struct named `NAME` according to specifications listed in the `RECIPE`.

## Customisation of a struct
The `RECIPE` may contain information pertinent to five distinct categories:
- Specification of `:dimensions`
- Specification of `:parameters` (their names and datatypes)
- Specification of `:auxiliary_info`
- Specification of `:conjugate` updates
- `:additional` information
Each type needs to be announced to `julia` by starting the list with the corresponding `Symbol` (or `QuoteNode`).
### Example
It's best to look at a simple example. Consider the definition of a Lorenz system:
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
The macro above expands to:
```julia
struct Lorenz{T} <: DiffusionProcess{Float64,3,3,UnboundedStateSpace()}
    p1::T
    p2::T
    p3::T
    σ::Float64
    function Lorenz(p1::T, p2::T, p3::T, σ::Float64) where T
        new{T}(p1, p2, p3, σ)
    end
end
```
which defines a [parametric type](https://docs.julialang.org/en/v1/manual/types/#Parametric-Types-1) `Lorenz{T}`, together with some handy auxiliary functions specific to any instance of `Lorenz{T}`. We may now instantiate the newly defined struct as in
```julia
P_f64 = Lorenz(1.0, 2.0, 3.0, 0.2)
P_SVectors = Lorenz([rand(DD.ℝ{3}) for i in 1:3]..., 0.2)
```
We can also call some functions auto-generated for the newly defined `Lorenz` struct, for instance
```julia
DD.param_names(P_f64) == (:p1, :p2, :p3, :σ)
DD.parameters(P_f64) == (1.0, 2.0, 3.0, 4.0)
```
There are in fact more functions that are being defined, see [Convenience Functions](/overview/convenience_functions/) for all functions that are being defined for each generated struct.

### More systematic explanations
The following information can be specified when using `@diffusion_process`.
#### `:dimensions`
Specification of the dimension of the process and the dimension of the driving Brownian motion. Must be written in a format: `process --> dimension` OR `wiener --> dimension` (eg. `process --> 4`)
#### `:parameters`
List of all parameters that the law depends on. Must be in one of the following formats (anywhere, if the parameter-name is set to `_` then uses default names based on a letter `p`):
- `single-parameter-name --> single-data-type` (eg. `σ --> Float64` or `σ --> T` (if `T` is a template argument) or `_ --> Vector{Int32}`). This defines a single parameter.
- `single-parameter-name --> (multiple-data-types,)` (eg. `σ --> (Float64, Int64)`, `_ --> (Int32, T)`). This defines as many parameters as there are specified types and appends the names with numbers to disambiguate multiple parameters with the same names.
- `single-parameter-name --> (number-of-parameters, data-type)` (eg. `σ --> (3, Float64)`) defines `number-of-parameters`-many parameters of the same data type.
- `(multiple-parameter-names,) --> single-data-type` (eg. `(α, β, γ) --> Int64`) defines as many parameters as there are names specified and sets them to be of the same type.
- `(multiple-parameter-names,) --> (multiple-data-types,)` (eg. `(α, β, γ) --> (Float64, Int64, T)`) defines as many parameters as there are names specified (there should be an equal number of names as there are types) and sets them to be of corresponding type.
#### `:auxiliary_info`
Information about the end-points of the diffusion. This is a useful feature for guided proposals or simulation of diffusion bridges, where the process is conditioned to hit certain end-point. The following fields can be defined:
- `t0`: the starting time-point
- `T`: the final time-point
- `x0`: the starting position
- `vT`: the observation at the terminal time
- `xT`: the state of the process at the terminal time
Each one of these fields can be defined in a format `field-name --> field-type`, for instance `T --> Float64`. **IMPORTANT** If at least one of the fields above is defined, then the field `xT` will be defined automatically and the diffusion constructor will default `xT` to a zero vector (or SVector or a Number, it will make a reasonable guess, but it may be wrong). The reason for this is to make employment of a blocking technique in other packages a bit easier.
#### `:conjugate`
This is a section need to be present if conjugate Gaussian updates are to be made in the setting of Bayesian inference. Three pieces of information can be specified:
- Function `phi`
- Function `nonhypo` in a format `nonhypo(x) --> non-smooth-coordinates-of-x` which specifies which coordinates of `x` have non-degenerate noise on them.
- `num_non_hypo` which specifies how many coordinates are with non-degenerate noise on them [TODO change to a simple count of `nonhypo` output length]
#### `:additional`
The additional information provides some additional decorators that helps the compiler use specialised functions when called on instances of corresponding diffusion processes. The following paramters can be specified
- `constdiff --> true` (or `false`) depending on whether the volatility coefficient is independent from the state variable (`false` by default)
- `linear --> true` (or `false`) to indicate that a diffusion has a linear structure (`false` by default).
- `statespace --> type-of-state-space-restrictions` (eg. `statespace --> UnboundedStateSpace()`) indicates any restrictions made on the state-space of a diffusion process.
- `eltype --> type-of-parameter` (eg. `eltype --> Float64`) disambiguate the parameter type in case multiple types are used. This is useful for automatic differentiation where the derivatives of only a subset of parameters are taken and it is the eltype of those parameters that is of interest.
