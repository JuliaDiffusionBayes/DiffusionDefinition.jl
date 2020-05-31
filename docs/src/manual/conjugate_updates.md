# Conjugate updates
****
In an MCMC setting it is sometimes possible to employ certain types of efficient updates called "conjugate updates". They can be employed if the likelihood and the prior fit together like a glove, allowing to pin down the posterior density by identifying it with one of the standard distributions. For the parameters appearing in the diffusion's drift such prior–posterior pair can **sometimes** be taken to be Gaussian–Gaussian.

This package does not deal with MCMC sampling for diffusions; however, it defines macros that automatically build functions that are required for conjugate updates in [DiffusionMCMC.jl](https://github.com/JuliaDiffusionBayes/DiffusionMCMC.jl). Below, we explain only how the macros introduced in DiffusionDefinition.jl work and give no further details about conjugate updates. To see more details about on this topic see [...].

## Conjugate Gausian updates
------
To perform conjugate Gaussian updates in [DiffusionMCMC.jl](https://github.com/JuliaDiffusionBayes/DiffusionMCMC.jl) the following functions need to be defines for your diffusion law:
```@docs
DiffusionDefinition.phi
DiffusionDefinition.num_non_hypo
DiffusionDefinition.hypo_a_inv
DiffusionDefinition.nonhypo
DiffusionDefinition.ignore_for_cu
```


In this package we implement a macro that makes this definition process more convenient
```@docs
DiffusionDefinition.@conjugate_gaussian
```

### Example

For instance, for the FitzHugh–Nagumo model in the `conjugate`-parameterization form from [here](@ref conjugate_fitzhugh_nagumo), we can define functions for the conjugate updates as follows:

```julia
@conjugate_gaussian FitzHughNagumoConjug begin
    :intercept --> (-x[2],)
    :ϵ --> (x[1]-x[1]^3+(1-3*x[1]^2)*x[2],)
    :s --> (one(x[1]),)
    :γ --> (-x[1],)
    :β --> (-one(x[1]),)
    :hypo_a_inv --> 1.0/P.σ^2
    :nonhypo --> 2:2
end
```
