<h1 align="center">
  <br>
  <a href="https://juliadiffusionbayes.github.io/DiffusionDefinition.jl/dev/"><img src="https://raw.githubusercontent.com/JuliaDiffusionBayes/DiffusionDefinition.jl/master/docs/src/assets/logo.png" alt="DiffusionDefinition.jl" width="200"></a>
  <br>
  DiffusionDefinition.jl
  <br>
</h1>

> A collection of convenient methods for defining diffusion processes and sampling from their laws.

<p align="center">
  <a href="https://JuliaDiffusionBayes.github.io/DiffusionDefinition.jl/stable">
    <img src="https://img.shields.io/badge/docs-stable-blue.svg"
         alt="Stable">
  </a>
  <a href="https://JuliaDiffusionBayes.github.io/DiffusionDefinition.jl/dev"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Dev"></a>
  <a href="https://travis-ci.com/JuliaDiffusionBayes/DiffusionDefinition.jl">
      <img src="https://travis-ci.com/JuliaDiffusionBayes/DiffusionDefinition.jl.svg?branch=master" alt="Build Status">
  </a>
</p>

<p align="center">
  <a href="#key-features">Key Features</a> •
  <a href="#installation">Installation</a> •
  <a href="#how-to-use">How To Use</a> •
  <a href="#related">Related</a> •
  <a href="#license">License</a>
</p>

## Key Features

- Convenient methods facilitating
  - defining diffusion laws
  - forward-sampling their trajectories
  - computing functionals of sampled paths
  - computing gradients of functionals of sampled paths with respect to diffusion parameters or with respect to the starting point of the trajectory
- A number of predefined diffusion processes that can be immediately loaded in and experimented on

## Installation

The package is not yet registered. To install it, type in:
```julia
] add https://github.com/JuliaDiffusionBayes/DiffusionDefinition.jl
```

## How To Use

See [the documentation](https://juliadiffusionbayes.github.io/DiffusionDefinition.jl/dev/).

## Related

DiffusionDefinition.jl belongs to a larger suite of packages in [JuliaDiffusionBayes](https://github.com/JuliaDiffusionBayes) designed to facilitate Bayesian inference for diffusion processes. Other packages in this suite consist of:
- [ObservationSchemes.jl](https://github.com/JuliaDiffusionBayes/ObservationSchemes.jl): a systematic way of encoding discrete-time observations for stochastic processes
- [GuidedProposals.jl](https://github.com/JuliaDiffusionBayes/GuidedProposals.jl): defining and sampling conditioned diffusion processes
- [ExtensibleMCMC.jl](https://github.com/JuliaDiffusionBayes/ExtensibleMCMC.jl): a modular implementation of the Markov chain Monte Carlo (MCMC) algorithms
- [DiffusionMCMC.jl](https://github.com/JuliaDiffusionBayes/DiffusionMCMC.jl): Markov chain Monte Carlo (MCMC) algorithms for doing inference for diffusion processes

## License

MIT
