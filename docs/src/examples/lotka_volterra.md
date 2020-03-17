# Lotka-Volterra
A simple, scalar-valued predator-prey model.

```math
\begin{align*}
\dd X_t &= (\alpha X_t - \beta X_t Y_t) \dd t + \sigma_1 \dd W^{(1)}_t \\
\dd Y_t &= (\delta X_t Y_t -\gamma Y_t)\dd t + \sigma_2 \dd W^{(2)}_t.
\end{align*}
```
Can be called with
```julia
@load_diffusion :LotkaVolterra
```

### Auxiliary diffusion
We also provide a linear diffusion that is obtained from linearising SDE above at the equilibrium point. This process can be used as an auxiliary diffusion in the setting of **guided proposals**. It solves the following stochastic differential equation
```math
\begin{align*}
\dd \widetilde{X}_t &= \left(- \frac{\beta\gamma}{\delta} \widetilde{Y}_t + \frac{\gamma\alpha}{\delta}\right) \dd t + \sigma_1 \dd W^{(1)}_t \\
\dd \widetilde{Y}_t &= \left(\frac{\alpha\delta}{\beta} \widetilde{X}_t-\frac{\alpha\gamma}{\beta}\right)\dd t + \sigma_2 \dd W^{(2)}_t,
\end{align*}
```
and can be called with
```julia
@load_diffusion :LotkaVolterraAux
```
