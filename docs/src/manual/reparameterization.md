# Reparameterizations
*****
Instances of `DiffusionProcess` are **mutable** so their fields may be changed directly. We provide additional convenience functions for reparameterizations. The first one clones the object and creates a new one with new parameters:
```@docs
DiffusionDefinition.clone
```
However, it's not recommended to rely on it and it will become deprecated in the near future.

A more efficient function, suitable for the MCMC setting is:
```@docs
DiffusionDefinition.set_parameters!
```
