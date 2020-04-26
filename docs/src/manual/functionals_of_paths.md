# Computing path-functionals when sampling
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
Y_0^{[3]}e^{-\int_0^T (Y_t^{[1]})^2e^{-t}\dd t + Y_T^{[2]}/100 }
```
for `Y` a path drawn under the Lorenz law, and $Y^{[i]}$ representing its `i`th coordinate. We may define functions `foo` as follows:
```julia
foo(P::Lorenz, y) = log(y[3])
foo(f_accum, P::Lorenz, y, t, dt, dW, i) = f_accum - exp(t)*y[1]^2*dt
foo(f_accum, P::Lorenz, y, ::Val{:final}) = exp(f_accum + y[2]/100.0)
X, foo_result = rand(P, tt, y1; f=foo)
```
## Gradients of path-functionals when sampling
Why would we care about computing path-functionals at the time of sampling instead of simply doing so after the path has been sampled and stored in a `Trajectory`? The answer is: we can leverage Julia's automatic differentiation implemented in a package [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) to compute gradients (and Hessians) of said functionals with respect to, say, a starting point or diffusion's parameters in a very efficient manner and without having to strain the memory by expanding the DataType used in `Trajectory`'ies.

!!! warning
    The syntax below might be simplified in the future to make it more user-friendly

### Example
For instance, in the example above, if we wanted to compute gradients with respect to diffusion parameters we can do the following:
```julia
# need a wrapper around foo
function foobar(θ)
    P = Lorenz(θ...)
    _, foo_result = DD.solve!(X, W, P, y1; f=foo)
    foo_result
end
tt = 0.0:0.01:1.0
X, W = trajectory(P, tt)
rand!(Wiener(), W)
y1 = SVector{3,Dual{Tag{typeof(foobar),Float64},Float64,4}}(
    Dual{Tag{typeof(foobar),Float64}}(rand(),0.0,0.0,0.0,0.0),
    Dual{Tag{typeof(foobar),Float64}}(rand(),0.0,0.0,0.0,0.0),
    Dual{Tag{typeof(foobar),Float64}}(rand(),0.0,0.0,0.0,0.0),
) # this is a little user-unfriendly
grad_at_θ° = ForwardDiff.gradient(foobar, θ°)
```
The trajectories `X` and `W` are of regular, default datatype specified by `P` (and do not use `Dual`s from `ForwardDiff`). Currently, `y1` needs to be already of the type that uses `Dual` numbers (that's how `solve!` knows what DataType to use) but this will most likely be simplified in the future.
