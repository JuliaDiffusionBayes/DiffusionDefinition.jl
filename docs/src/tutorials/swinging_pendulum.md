# [Defining a diffusion process for a swinging pendulum](@id tutorials_start)
***
> In this tutorial we will define a diffusion and simulate trajectories from it.

## The model
----
We start with the differential equation for the angular position of a swinging pendulum, which is given by the second
order ordinary differential equation
```math
 \frac{\dd^2 x(t)}{\dd t^2} + θ^2 \sin(x(t)) = 0.
```
Here $x(t)$ gives the angular position at time $t$ and $θ$ is  the angular velocity of the linearised pendulum.
Under the assumption that the acceleration is in fact a white-noise process, we obtain the
Stochastic Differential Equation (SDE)
```math
 \begin{align}
 \dd X_t &=\begin{bmatrix} 0 & 1 \\ 0 & 0 \end{bmatrix} X_t \dd
 t + \begin{bmatrix} 0 \\ -θ^2 \sin(X_{1t})\end{bmatrix} \dd t +
  \begin{bmatrix} 0 \\ \gamma \end{bmatrix} \dd W_t, \end{align}
```
where $X_t=\begin{bmatrix} X_{t1} & X_{t2}\end{bmatrix}^\prime$.


## Defining diffusion
----
To define the diffusion for the pendulum, we use the `DiffusionDefinition.jl`-package.
```julia
using Plots, DiffusionDefinition, StaticArrays
```
To avoid typing `DiffusionDefinition` fully, we define
```julia
const DD = DiffusionDefinition;
```
To define the diffusion, we specify the dimensions of the state space and Wiener noise and the parameters appearing in the SDE using the `diffusion_process`-macro. In the call to this macro we specify the dimensions of the state space of the diffusion and driving Wiener process, as also parameters appearing in their definition.
```julia
@diffusion_process Pendulum begin
    :dimensions
    process --> 2
    wiener --> 1

    :parameters
    (θ, γ) --> (2, Float64)
end
```
Note the printed message and that by issuing the command `?Pendulum` information on the constructed `struct` is displayed.
Suppose $θ=2$ and $γ=0.5$, then we can define an instance of the `struct` called `Pendulum` by setting
```julia
P = Pendulum(2.0, 0.5)
```
There are various convenience functions, among which
```julia
DD.parameters(P)
DD.parameter_names(P)
```
The general form in which the diffusion is specified uses the notation
```math
 \dd X_t = b(t,X_t)\dd t + \sigma(t,X_t) \dd W_t
```
So the drift is denoted by $b$ and the diffusivity by $σ$ (this is the notation in the well known books by Rogers and Williams on stochastic calculus).
As the functions $b$ and $σ$ are not exported in the `DiffusionDefinition`-package (to avoid conflicts with other packages), these need to
be defined using `DiffusionDefinition.b` and `DiffusionDefinition.σ`.
```julia
DD.b(t, x, P::Pendulum) = @SVector [x[2],-P.θ^2 * sin(x[1])]
DD.σ(t, x, P::Pendulum) = @SMatrix [0.0 ; P.γ]
```
While not strictly necessary, it is preferable to use static vectors and matrices using commands from the `StaticArrays`-package.

## Sampling trajectories
----
Now we can sample a trajectory of the diffusion. For that, we need a starting point and grid (on which the solution to the SDE is approximated using Euler discretisation).
```julia
tt = collect(0.0:0.005:8.0)  # grid
x0 = @SVector [1.0, 0.0]  # starting point
X = DD.rand(P, tt, x0);
```
Now `X` is a trajectory, with `X.t` denoting the timegrid, and `X.x` the simulated trajectory.

# Plotting
----
Next, we plot the simulated trajectories
```julia
#using PyPlot
pyplot()
plot(X, Val(:vs_time))
```
