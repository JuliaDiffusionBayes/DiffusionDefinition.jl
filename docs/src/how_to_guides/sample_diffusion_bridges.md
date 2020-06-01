# How to sample diffusion bridges?
***
A diffusion bridge is a diffusion conditioned on its end-point. A diffusion bridge is a special type of a conditioned diffusion where the conditioning is on its end point.

The most convenient, efficient and robust method of sampling conditioned diffusions (including diffusion bridges) is to use the package [GuidedProposals.jl](https://juliadiffusionbayes.github.io/GuidedProposals.jl/dev/) that has been designed to address precisely this problem.

Alternatively, you may follow the [tutorial on sampling diffusion bridges using  rejection and importance sampling with Wiener proposals](@ref tutorial_sampling_diff_bridges). In vast majority of cases however, this last method will not be suitable for production models.
