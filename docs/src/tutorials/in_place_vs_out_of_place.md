# Understanding in-place and out-of-place methods and how to call them
****
> [DiffusionDefinition.jl](https://juliadiffusionbayes.github.io/DiffusionDefinition.jl/dev/) provides support for both in-place and out-of-place sampling. So which one should you use? In this tutorial we will illustrate which one of the two modes should be used on the example of a Lorenz96 system.

## Lorenz '96 system
---
In this tutorial we will be using a Lorenz '96 system. More details about it can be found [here](@ref lorenz96_system). We will be using it because it is a model that can be defined for different dimension $D$ of the state space. The model is pre-defined in the package and can be simply loaded in with `@load_variable_diffusion`, however, for didactic purposes we will define it explicitly using `@diffusion_process`.

## Out-of-place methods
---
Out-of-place methods are the type of methods that are allowed to allocate new space for performing computations. If you are able to ascertain that your variables live on a `stack allocated memory`, then out-of-place operations can be very efficient. For [DiffusionDefinition.jl](https://juliadiffusionbayes.github.io/DiffusionDefinition.jl/dev/) this corresponds to the state of the diffusion process being represented either by `Number`s or `SArray`s. So the question becomes: if the state-space of your diffusion is $\RR^d$, when does it make sense to have it represented by `SArray`s? Should we simply always use `SArray`s? To answer this question you should simply consult the documentation of [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl). It states clearly that if your objects become greater in dimension than $100$, then the benefits of using `SArray`s start to taper off and it is better to use regular arrays.

For diffusion processes the largest element will most often be a volatility matrix. If it is represented by a dense array, then for diffusions living on $\RR^{10}$ and higher we should start re-considering use of out-of-place methods. If the volatility matrix can be kept sparse (or in particular diagonal) then the limit for using out-of-place methods are rather diffusions living in $\RR^{100}$. As a general rule of thumb you should use `SArray`s (or `Number`s) for dimension less than that and regular `Array`s otherwise.


### Example

Let's define a 50 dimensional Lorenz '96 system:
```julia
@diffusion_process Lorenz96_OOP begin
    :dimensions
    process --> 50
    wiener --> 50

    :parameters
    (θ, σ) --> Float64

    :additional
    constdiff --> true
end

const si1 = SVector{50,Int64}(mod1.(2:51, 50))
const si2 = SVector{50,Int64}(mod1.(-1:48, 50))
const si3 = SVector{50,Int64}(mod1.(0:49, 50))
const si4 = SVector{50,Int64}(1:50)

@inline function b(t, x, P::Lorenz96_OOP)
    ( view(x, si1) .- view(x, si2) ) .* view(x, si3) .- view(x, si4) .+ P.θ
end

σ(t, x, P::Lorenz96_OOP) = P.σ * SDiagonal{50,Float64}(I)
```
Note that any subindexing must be done with `SVector`s to avoid creating regular arrays first that are later converted to `SArray`s.

## In-place methods
---
In-place methods are those assumed to have a pre-allocated, fixed amount of memory on which all operations are to be performed. Based on the previous section, you can imagine that in-place methods should simply be used whenever out-of-place methods cannot.


### Example

For instance if we bump the dimension of the Lorenz '96 system to, say `10000`, then we should most certainly switch to regular arrays. This time however, the output of the drift and volatility should be saved into a pre-allocated `buffer`.

```julia
@diffusion_process Lorenz96_IP begin
    :dimensions
    process --> 10000
    wiener --> 10000

    :parameters
    (θ, σ) --> Float64

    :additional
    constdiff --> true
end

const i1 = mod1.(2:10001, 10000)
const i2 = mod1.(-1:9998, 10000)
const i3 = mod1.(0:9999, 10000)
const i4 = 1:10000

function b!(buffer, t, x, P::Lorenz96_IP)
    buffer.b .= ( view(x, i1) .- view(x, i2) ) .* view(x, i3) .- view(x, i4) .+ P.θ
end

σ!(buffer, t, x, P::Lorenz96_IP) = (buffer.σ.diag .= P.σ)
```
