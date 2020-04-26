# Sine diffusion
A simple, scalar diffusion process displaying multimodality on a path space. It solves the following stochastic differential equation
```math
\dd X_t = (a + b\sin(c X_t))\dd t + \sigma \dd W_t,\quad t\in[0,T],\quad X_0=x_0.
```
It can be called with
```julia
@load_diffusion :Sine
```

### Auxiliary diffusion
We define an additional, linear diffusion that can be used in the setting of **guided proposals**. It solves the following SDE
```math
\dd \widetilde{X}_t = \left(\frac{x_T-x_0}{T} + \frac{t}{5T}\widetilde{X}_t\right)\dd t + \sigma \dd W_t,\quad t\in[0,T],\quad X_0=x_0,
```
and can be called with
```julia
@load_diffusion :SineAux
```
