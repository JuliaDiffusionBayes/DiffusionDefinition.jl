# How to compute expected values of gradients of functionals?
***
To compute expectations of gradients of functionals we may simply use Monte Carlo simulations. For instance, suppose we have a Lorenz system
```julia
@load_diffusion Lorenz
θ° = [10.0, 28.0, 8.0/3.0, 1.0]
tt = 0.0:0.001:10.0
y1 = @SVector [-10.0, -10.0, 25.0]
```

and the following functional to compute
```julia
foo(P::Lorenz, y) = 0.0
foo(f_accum, P::Lorenz, y, t, dt, dW, i) = f_accum - exp(-t)*y[1]^2/1000.0*dt
foo(f_accum, P::Lorenz, y, ::Val{:final}) = exp(f_accum)
```

Then, we can use, say, function `grad_θ` that computes gradients with respect to the parameters and embed it in a Monte Carlo scheme as follows:

```julia
function expectation_of_grad(θ°, y1, W, X, Law, f, num_rep=1e4)
    grads = map(1:num_rep) do i
        rand!(Wiener(), W)
        # gradient computation
        DD.grad_θ(θ°, y1, W, X, Law, f)
    end
    E = sum(grads)/num_rep
    C = mapreduce(g->(g-E)*(g-E)', +, grads)/num_rep
    E, C
end
```

We can now compute the expected values of gradients together with the covariance:

```julia
X, W = trajectory(tt, Lorenz(θ°...))
E_grad, C_grad = expectation_of_grad(θ°, y1, W, X, Lorenz, foo)
```

which for this particular case yields:
```julia
julia> E_grad
4-element SArray{Tuple{4},Float64,1,4} with indices SOneTo(4):
  9.43393204196648e-5
 -0.0029595272871572606
 -0.024297701677402275
  0.0001291653982423314

julia> C_grad
4×4 SArray{Tuple{4,4},Float64,2,16} with indices SOneTo(4)×SOneTo(4):
  4.99701e-6   1.0669e-5    -3.05976e-5    4.25573e-5
  1.0669e-5    4.34692e-5   -0.000240906   0.000171052
 -3.05976e-5  -0.000240906   0.00191821   -0.0010401
  4.25573e-5   0.000171052  -0.0010401     0.000800068
```
