# Lorenz '96 system
A model of atmospheric convection where each coordinate of the $D$-dimensional process $X$ corresponds to a position on a periodic lattice that is supposed to be a proxy for a latitude circle on Earth. The diffusion $X$ is defined as a solution to the following system of SDEs:
```math
\dd X_t = \left[\left(X^{(i+1)}_t-X^{(i-2)}_t\right)X^{(i-1)}_t-X^{(i)}_t+F\right]\dd t + \sigma_i \dd W^{(i)}_t,\quad t\in[0,T],\quad X^{(i)}_0=x^{(i)}_0,
```
with $i\in\{1,\dots,D\}$ a cycling  index.
The process can be called with
```julia
@load_variable_diffusion :Lorenz96 D
```
where $D$ is a positive integer, indicating chosen dimension.
