# Sampling diffusion trajectories
This package extends the functionality of `Random.rand!` and `Base.rand` to sampling of trajectories of diffusion processes. This package does not deal with the most versatile implementation of the unconditioned SDE samplers (for that see [DifferentialEquations.jl](https://docs.sciml.ai/stable/)) or conditioned SDE samplers (for that see [Bridge.jl](https://github.com/mschauer/Bridge.jl)). Here we provide only the necessary (and efficient) functionality required for imputing diffusion processes and conducting Bayesian inference for them.

The containers for sampled trajectories are chosen to be instances of `Trajectory` from the package [Trajectories.jl](https://github.com/mschauer/Trajectories.jl). The functions exported by `Trajectories.jl` are re-exported by package `DiffusionDefinition.jl`.

### Copied from readme:
For instance, to sample a three-dimensional standard Brownian motion use:
```julia
const DD = DiffusionDefinition
tt = collect(0.0:0.01:1.0)
wiener_path = rand(Wiener(), zero(DD.ℝ{3}), tt)
```
The wiener path can then be used in an Euler-Maruyama scheme to compute a trajectory under a given diffusion law:
```julia
XX = trajectory(tt, DD.ℝ{3})
P = Lorenz(28.0, 10.0, 8.0/3.0, 2.0)
DD.solve!(XX, wiener_path, P, zero(DD.ℝ{3}))
```
The two-step sampling recipe above can also be done at once by calling:
```julia
XX = rand(P, zero(DD.ℝ{3}), tt)
```
