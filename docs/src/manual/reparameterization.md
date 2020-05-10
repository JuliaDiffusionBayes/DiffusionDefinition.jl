# Reparameterizations
Instances of `DiffusionProcess` are immutable, and usually, so are their fields with parameters. Consequently, reparameterizations almost always involve constructing new objects.

In an MCMC setting reparameterizations of diffusion laws happen very frequently and thus we provide a convenience function suitable for that setting:
```@docs
DiffusionDefinition.clone
```

[TODO make it less convoluted]
