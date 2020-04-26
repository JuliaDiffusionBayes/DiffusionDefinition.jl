# DiffusionDefinition.jl

This is a compact utility package for defining diffusion processes and sampling from their laws. It is created to work in conjunction with a suite of packages in [JuliaDiffusionBayes](https://github.com/JuliaDiffusionBayes) that provide tools for Bayesian inference for diffusion processes. However, it can also be used on its own to:
- define diffusion laws
- forward-sample their trajectories
- compute functionals of sampled paths
- compute gradients of functionals of sampled paths with respect to diffusion parameters or with respect to the starting point of the trajectory

Depending on what is your level of experience and, most importantly, your intended use of this package you might consider from starting at different places:

- For a quick overview of this package's main functionality see [Get started](@ref get_started).
- For a systematic introduction to all functionality introduced in this package see the [Manual](@ref defining_diffusion)
- For a didactic introduction to problems that can be solved using our package see the [Tutorials](@ref tutorials_start)
- If you have a problem that you think can be addressed with this package, then check out the [How-to guides](@ref how_to_guides) to see if the answer is already there.
