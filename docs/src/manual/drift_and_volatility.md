# Defining the drift and the volatility coefficient
To complete the definition of a diffusion process we need to specify its drift, as well as its volatility coefficient. For the `Lorenz` example from the [previous section](@ref Lorenz_example) we can do this by writing:
```julia
const DD = DiffusionDefinition
using StaticArrays

function DD.b(t, x, P::Lorenz)
    @SVector [
        P.p1*(x[2]-x[1]),
        P.p2*x[1] - x[2] - x[1]*x[3],
        x[1]*x[2] - P.p3*x[3]
    ]
end

DD.σ(t, x, P::Lorenz) = SDiagonal(P.σ, P.σ, P.σ)
```
The two functions above are **out-of-place** i.e. they return new vectors (that live on a stack, because of `StaticArrays`). We may alternatively define the drift and diffusion coefficients to be **in-place** as follows:
```julia
function DD.b!(buffer, t, x, P::Lorenz)
    buffer.b[1] = P.p1*(x[2]-x[1])
    buffer.b[2] = P.p2*x[1] - x[2] - x[1]*x[3]
    buffer.b[3] = x[1]*x[2] - P.p3*x[3]
end

DD.σ!(buffer, t, x, P::Lorenz) = (buffer.σ.diag .= P.σ)
```
In this case case the output is saved to a `buffer`, which must have appropriate fields `b` and `σ` with enough pre-allocated space (see also [the section on buffers](@ref explain_buffers)).
!!! note
    * all of the functions `DD.b`, `DD.b!`, `DD.σ` and `DD.σ!` are defined to overload the functionality inside the `DiffusionDefinition` module (accessing it via `DD`) and **NOT** the `Main` module.
    * the arguments for the out-of-place method are `(t, x, P::DIFFUSION_NAME)`, whereas those for in-place are `(buffer, t, x, P::DIFFUSION_NAME)`.

!!! tip
    Always use `StaticArrays` for out-of-place drift and volatility! If functions using `DD.b` and `DD.σ` are faster with regular arrays, then you shouldn't be using out-of-place methods in the first place, but `DD.b!` and `DD.σ!` instead. A general rule of thumb is to use `DD.b` and `DD.σ` for low dimensional diffusions (up to dimension `~10` for elliptic diffusions with dense volatility coefficients or up to dimension `~100` for those with sparse volatility coefficients) and use in-place methods otherwise.

## [Telling Julia whether to use (`DD.b`, `DD.σ`) or (`DD.b!`, `DD.σ!`)](@id default_types_for_P)
Some functions implemented in this package have two versions: one relying on out-of-place methods, another on in-place methods. For instance:
```@docs
DiffusionDefinition.solve!(XX, WW, P, y1)
```
and
```@docs
DiffusionDefinition.solve!(XX, WW, P, y1, buffer)
```
Calling one or the other will tell Julia which pair of functions (`DD.b`, `DD.σ`) or (`DD.b!`, `DD.σ!`) to use. However, most functions do not explicitly come in two versions, and instead, they rely on some hint from the user to decide on the mode of computation. There are two ways in which Julia can be given such hints:
- If the function accepts optional inputs (such as a starting point) or expects to receive containers that the results are saved into (such as containers for diffusion paths), then the `DataType` of said inputs will be used to decide on the mode of computation (if `DataType` used for state space is immutable, then out-of-place methods are used, otherwise in-place methods are used). These hints will overwrite the second method below.
- Sometimes no such inputs can be passed or there is an option of relying on defaults. In that case `Julia` will use default information specified by functions:
```@docs
DiffusionDefinition.default_type
DiffusionDefinition.default_wiener_type
```
By default, these two are set to `StaticArrays` resulting in out-of-place computations by default. To change that, overwrite the two functions for your diffusion type. For instance, for the `Lorenz` example, to change the default mode of computation to out-of-place write:
```julia
default_type(::Lorenz) = Vector{Float64}
default_wiener_type(::Lorenz) = Vector{Float64}
```
