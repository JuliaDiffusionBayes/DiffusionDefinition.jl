# Lorenz '63 system
Famous Lorenz attractor, a three-dimensional elliptic diffusion, a solution to the following stochastic differential equation
```math
\begin{align*}
\dd X_t &= \theta_1 (Y_t - X_t) \dd t + \sigma \dd W^{(1)}_t \\
\dd Y_t &= [X_t (\theta_2 - Z_t) - Y_t]\dd t + \sigma \dd W^{(2)}_t\\
\dd Z_t &= [X_t Y_t - \theta_3 Z_t]\dd t + \sigma \dd W^{(3)}_t.
\end{align*}
```
Can be imported with the following command
```julia
@load_diffusion :Lorenz
```
### Auxiliary diffusion
We additionally provide an implementation of a linear diffusion that can be used in a setting of **guided proposals**. It is defined as a solution to the following SDE:
```math
\begin{align}
\dd \wt{X}_t &= \left[-\theta_1 \wt{X}_t + \theta_2\wt{Y}_t\right]\dd t + \sigma \dd W^{(1)}_t,\\
\dd \wt{Y}_t &= \left[ (\theta_2-z_T)\wt{X}_t - \wt{Y}_t - x_T\wt{Z}_t + x_Tz_T \right]\dd t + \sigma \dd W^{(2)}_t,\\
\dd \wt{Z}_t &= \left[ y_T\wt{X}_t x_T\wt{X}_t -\theta_3\wt{Z}_t -x_Ty_T \right]\dd t + \sigma \dd W^{(3)}_t.
\end{align}
```
It can be called with
```julia
@load_diffusion :LorenzAux
```
