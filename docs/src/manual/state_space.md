# [Restricting diffusion's state space](@id state_space_restrictions)
******
Many standard diffusions are defined only on a subset of $$\RR^d$$, and thus, it is often needed to restrict the state space on which trajectories can be sampled. We provide a set of structs that may be passed at the time of defining a diffusion law (calling the macro `@diffusion_process`) to restrict the state space of the diffusion
```@docs
DiffusionDefinition.UnboundedStateSpace
DiffusionDefinition.LowerBoundedStateSpace
DiffusionDefinition.UpperBoundedStateSpace
DiffusionDefinition.BoundedStateSpace
```
