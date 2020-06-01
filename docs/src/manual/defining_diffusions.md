# [Defining diffusion](@id defining_diffusion)
--------------------------------------------------------------------------------
The main utility macro introduced in this package is `@diffusion_process`. It facilitates very concise definitions of structs characterizing diffusion processes.

There are two (optionally three) parts expected by the macro:
- the name of the diffusion (which, optionally, may contain template parameters in the curly brackets)
- and the recipe for defining a struct
```julia
@diffusion_process NAME{TEMPLATE_PARAMETERS} begin
  RECIPE
end
```
the snippet of code above creates a struct named `NAME` according to specifications listed in the `RECIPE`.

## Customisation of a struct
The `RECIPE` may contain information pertinent to five distinct categories:
- Specification of `:dimensions`
- Specification of `:parameters` (their names and datatypes)
- Specification of `:constant_parameters` (their names and datatypes)
- Specification of `:auxiliary_info`
- `:additional` information
Each type needs to be announced to `julia` by starting the list with the corresponding `Symbol` (or `QuoteNode`).

!!! note
    For many users knowing only about `:dimensions`, `:parameters` and `:additional` will be sufficient and the other categories will not be of much importance.

### [Example](@id Lorenz_example)
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
    diagonaldiff --> true
end
```
The macro above expands to:
```julia
struct Lorenz{T} <: DiffusionProcess{Float64, 3, 3, UnboundedStateSpace()}
    p1::T
    p2::T
    p3::T
    σ::Float64
    function Lorenz(p1::T, p2::T, p3::T, σ::Float64) where T
        new{T}(p1, p2, p3, σ)
    end
    function Lorenz(; p1::T, p2::T, p3::T, σ::Float64) where T
        new{T}(p1, p2, p3, σ)
    end
end
```
which defines a [parametric type](https://docs.julialang.org/en/v1/manual/types/#Parametric-Types-1) `Lorenz{T}`, together with some handy auxiliary functions specific to any instance of `Lorenz{T}`. We may now instantiate the newly defined struct as in
```julia
P_f64 = Lorenz(10.0, 28.0, 8.0/3.0, 0.2)
# or
P_f32 = Lorenz(10.0f0, 28.0f0, 8.0f0/3.0f0, 0.2)
```
We can also call some functions that were auto-generated for the newly defined `Lorenz` struct, for instance
```julia
DD.parameter_names(Lorenz) == (:p1, :p2, :p3, :σ)
DD.parameter_names(P_f64) == (:p1, :p2, :p3, :σ)
DD.parameters(P_f64) == Dict(:p1 => 10.0, :p2 => 28.0, :p3 => 8.0/3.0, :σ => 0.2)
```
More functions are automatically defined in the background for each generated `DiffusionProcess`, to learn more about them see the list of [Utility Functions](@ref utility_functions).

# Systematic explanations
--------------------------------------------------------------------------------
The following information can be specified in the definition of a diffusion law when calling a macro `@diffusion_process`.

!!! tip
    All keywords read by `@diffusion_process` are **case-insensitive**.

#### Category 1: `:dimensions`
Specification of the dimension of the process and the dimension of the driving Brownian motion. Must be written in a format:
- `process --> dimension` OR
- `wiener --> dimension` (eg.
  - `process --> 4`).
In both cases, if left unspecified then defaults to `dimension=1`.

!!! note "Alternative keywords"
    `:dimensions` is not the only keyword that will be recognized by `@diffusion_process` as declaring the dimensions of the process. Alternative names that could be used in place of `:dimensions` are: `:dim`, `:dims`, `:dimension`

#### Category 2: `:parameters`
List of all parameters that the law depends on.

!!! note
    `@diffusion_process` understands `_` as "*the user doesn't care about the name, so let's use a generic name based on the letter `p` and append it with a disambiguation number so that if there are more than one `p`'s they are not confused with each other*"

The parameters must be specified in one of the following formats:
- `single-parameter-name --> single-data-type`, (eg.
  - `σ --> Float64` or
  - `σ --> T` (if `T` is one of the template's labels) or
  - `_ --> Vector{Int32}`).
  This defines a single parameter.
- `single-parameter-name --> (multiple-data-types,)` (eg.
  - `σ --> (Float64, Int64)`,
  - `_ --> (Int32, T)`).
  This defines as many parameters as there are specified types and appends the names with numbers to disambiguate multiple parameters with the same names.
- `single-parameter-name --> (number-of-parameters, data-type)` (eg.
  - `σ --> (3, Float64)`)
  defines `number-of-parameters`-many parameters of the same data type.
- `(multiple-parameter-names,) --> single-data-type` (eg.
  - `(α, β, γ) --> Int64`)
  defines as many parameters as there are names specified and sets them to be of the same type.
- `(multiple-parameter-names,) --> (multiple-data-types,)` (eg.
  - `(α, β, γ) --> (Float64, Int64, T)`)
  defines as many parameters as there are names specified (there should be an equal number of names as there are types, it will throw an `AssertionError` otherwise) and sets them to be of the corresponding type.

!!! note "Alternative keywords"
    `:parameters` keyword alternatives: `:param`, `:params`

#### Category 3: `:constant_parameters`
These can be defined in exactly the same way as `:parameters`. The only purpose for splitting the parameters into `:constant_parameters` and `:parameters` is to indicate to `Julia` that the set of all parameters may be split into two conceptually different groups. In particular, `@diffusion_process` defines utility functions that act differently with `:constant_parameters` and `:parameters`:
```@docs
DiffusionDefinition.parameters
DiffusionDefinition.const_parameters
DiffusionDefinition.var_parameters
DiffusionDefinition.parameter_names
DiffusionDefinition.const_parameter_names
DiffusionDefinition.var_parameter_names
```
!!! note
    The split into `constant` and `variable` parameters is not done at a compile time and can also be done by hand after the `struct` with the diffusion has been constructed. For instance, in the `Lorenz` example above we have:

    ```julia
    julia> DD.parameters(P_f64)
    Dict{Symbol,Float64} with 4 entries:
      :p2 => 28.0
      :σ  => 0.2
      :p1 => 10.0
      :p3 => 2.66667

    julia> DD.const_parameters(P_f64)
    Dict{Any,Any} with 0 entries

    julia> DD.var_parameters(P_f64)
    Dict{Symbol,Float64} with 4 entries:
      :p2 => 28.0
      :σ  => 0.2
      :p1 => 10.0
      :p3 => 2.66667
    ```

    If we wanted to change our mind and define, say, `:p1` and `:p3` as `constant` parameters we could do that by overwriting the definition of the `const_parameter_names(::Type{<:CustomDiffusionLaw})` method, as all other functions of this type are computed as a byproduct of `parameter_names(::Type{<:CustomDiffusionLaw})` and `const_parameter_names(::Type{<:CustomDiffusionLaw})` and the former one should never change.

    ```julia
    DD.const_parameter_names(::Type{<:Lorenz}) = (:p1, :p3)
    ```

    That's all that needs to be changed.


!!! note "Alternative keywords"
    `:constant_parameters` keyword alternatives: `const_parameters`, `:const_param`, `:const_params`, `:constparameters`, `:constparam`, `:constparams`, `:constant_param`, `:constant_params`, `:constantparameters`, `:constantparam`, and `:constantparams`.

#### Category 4: `:auxiliary_info`
Information about the end-points of the diffusion. This is a useful feature for guided proposals or simulation of diffusion bridges, where the process is conditioned to hit a certain end-point. In particular, it is used quite extensively in the package [GuidedProposals.jl](https://github.com/JuliaDiffusionBayes/GuidedProposals.jl). The following fields can be defined:
- `t0` (also accepts `:t_0`): the starting time-point
- `T`: the final time-point
- `x0` (also accepts `:y0`, `:state0`, `:x_0`, `:y_0`, `:state_0`): the starting position
- `v0` (also accepts `:obs0`, `:v_0`, `:obs_0`): the starting observation
- `vT` (also accepts `:obsT`, `:v_T`, `:obs_T`): the observation at the terminal time
- `xT` (also accepts `:yT`, `:stateT`, `:x_T`, `:y_T`, `:state_T`): the state of the process at the terminal time
Each one of these fields can be defined in a format `field-name --> field-type` (e.g. `T --> Float64`).
!!! warning "IMPORTANT"
    If at least one of the fields above is defined, then the field `xT` will be defined automatically and the diffusion constructor will default `xT` to a zero `Vector` (or `SVector` or `Number`, it will make a reasonable guess, but it may be wrong). The reason for this is to make employment of a **blocking** technique in other packages a bit easier.

!!! note "Alternative keywords"
    `auxiliary_info` keyword alternatives: `:aux_info`, `:end_points`, `:end_point_info`

#### Category 5: `:additional`
The additional information provides some additional decorators that helps the compiler use specialized functions when called on instances of corresponding diffusion processes. The following information can be specified
- `constdiff --> true` (or `false`) depending on whether the volatility coefficient is independent from the state variable (`false` by default).
  - Alternative keywords: `:constvola`, `:constdiffusivity`, `:constvolatility`, `:constσ`, `:constantdiff`, `:constantvola`, `:constantdiffusivity`, `:constantvolatility`, `:constantσ`.
- `diagonaldiff --> true` (or `false`) to indicate that the volatility coefficient is represented by a diagonal matrix (`false` by default).
  - Alternative keywords: `:diagonalvola`, `:diagonaldiffusivity`, `:diagonalvolatility`, `:diagonalσ`, `:diagdiff`, `:diagvola`, `:diagdiffusivity`, `:diagvolatility`, `:diagσ`.
- `sparsediff --> true` (or `false`) to indicate that the volatility coefficient is a represented by a sparse matrix (`false` by default).
  - Alternative keywords: `:sparsevola`, `:sparsediffusivity`, `:sparsevolatility`, `:sparseσ`.
- `linear --> true` (or `false`) to indicate that a diffusion has a linear structure (`false` by default).
- `diagonalbmat --> true` (or `false`) to indicate that a `B` matrix of a linear diffusion (with a drift `b(x):=Bx+β`) is diagonal (`false` by default).
  - Alternative keywords: `:diagonalb`, `:diagonalbmatrix`.
- `sparsebmat --> true` (or `false`) to indicate that a `B` matrix of a linear diffusion (with a drift `b(x):=Bx+β`) is sparse (`false` by default).
  - Alternative keywords: `:sparseb`, `:sparsebmatrix`.
- `statespace --> type-of-state-space-restrictions` (eg. `statespace --> UnboundedStateSpace()`) indicates any restrictions made on the state-space of a diffusion process.
- `eltype --> type-of-parameter` (eg. `eltype --> Float64`) disambiguate the parameter type in case multiple types are used. This is useful for automatic differentiation where the derivatives of only a subset of parameters are taken and it is the eltype of those parameters that is of interest. [**TODO** come back, not sure anymore if it has any use]

!!! note "Alternative keywords"
    `:additional` keyword alternatives: `:extra`
