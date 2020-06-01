# Computing path-functionals when sampling
--------------------------------------------------------------------------------
`rand`, `rand!` and `solve!` accept an additional, named argument `f`—a functional of a path. By default it is turned off, but we may define some `Julia` function and pass it on to `rand`, `rand!` or `solve!` to be evaluated during the call to `solve!`. For out-of-place computations three versions of function need to be defined:
```julia
# called at the very start of solve!
f_accum = foo(P, y)
# called at the beginning of every iteration of the Euler-Maruyama scheme
f_accum = foo(f_accum, P, y, t, dt, dW, i)
# called at the very end, just before return statement
f_accum = foo(f_accum, P, y, Val(:final))
```
For in-place computations these three functions must have a slightly different form
```julia
# called at the very start of solve!
f_accum = foo(buffer, P, y) # the buffer needs to accommodate the needs of function f
# called at the beginning of every iteration of the Euler-Maruyama scheme
foo(buffer, P, y, tt[i-1], dt, i-1)
# called at the very end, just before return statement
foo(buffer, P, y, _FINAL)
```
If `foo` is appropriately defined, then we may pass it to, say, `rand` as follows:
```julia
X, foo_result = rand(P, tt, y1; f=foo)
```
Then, the output of `rand` is not only the sampled path `X`, but also the result of evaluating the functional defined by `foo`.
### Example
As an example, consider computing:
```math
Y_0^{[3]}e^{-\frac{1}{100}\int_0^T (Y_t^{[1]})^2e^{-t}\dd t + Y_T^{[2]}/100 }
```
for `Y` a path drawn under the Lorenz law, and $Y^{[i]}$ representing its `i`th coordinate. We may define functions `foo` as follows:
```julia
foo(P::Lorenz, y) = log(y[3])
foo(f_accum, P::Lorenz, y, t, dt, dW, i) = f_accum - exp(-t)*y[1]^2/100.0*dt
foo(f_accum, P::Lorenz, y, ::Val{:final}) = exp(f_accum + y[2]/100.0)
X, foo_result = rand(P, tt, y1; f=foo)
```
## Gradients of path-functionals when sampling
Why would we care about computing path-functionals at the time of sampling instead of simply doing so after the path has been sampled and stored in a `Trajectory`? The answer is: we can leverage Julia's automatic differentiation implemented in the package [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) to compute gradients (and Hessians) of said functionals with respect to, say, a starting point or diffusion's parameters in a very efficient manner and without having to strain the memory by expanding the DataType used in `Trajectory`'ies.

### Example
For instance, in the example above, if we wanted to compute gradients with respect to diffusion parameters we can do the following:
```julia
# need a wrapper around foo
function foobar(θ)
    P = Lorenz(θ...)
    _, foo_result = DD.solve!(X, W, P, y1_θ; f=foo)
    foo_result
end
# stating parameters
θ° = @SVector[10.0, 28.0, 8.0/3.0, 1.0]
tt = 0.0:0.001:10.0
y1 = @SVector [-10.0, -10.0, 25.0]
X, W = trajectory(tt, P)
rand!(Wiener(), W)
# gradient preparation
y1_θ_type = similar_type(θ°, Dual{Tag{typeof(foobar),eltype(y1)},eltype(y1),length(θ)}, Size(y1))
y1_θ = y1_θ_type(y1)

grad_at_θ° = ForwardDiff.gradient(foobar, θ°)
```
The trajectories `X` and `W` are of regular, default datatype specified by `P` (and do not use `Dual`s from `ForwardDiff`). `y1` needs to be already of the type that uses `Dual` numbers (that's how `solve!` knows what DataType to use).

The steps above are packaged in a convenience function
```@docs
DiffusionDefinition.grad_θ
```
so all the steps above are encapsulated by a call
```julia
# stating parameters
θ° = @SVector[10.0, 28.0, 8.0/3.0, 1.0]
tt = 0.0:0.001:10.0
y1 = @SVector [-10.0, -10.0, 25.0]
X, W = trajectory(P, tt)
rand!(Wiener(), W)

# gradient computation
grad_at_θ° = grad_θ(θ°, y1, W, X, Lorenz, foo)
```

Additionally, we provide:
```@docs
DiffusionDefinition.grad_y1
```
so that gradients with respect to a starting position may be computed as follows:
```julia
P = Lorenz(θ°...)
DD.grad_y1(y1, W, X, P, foo)
```
