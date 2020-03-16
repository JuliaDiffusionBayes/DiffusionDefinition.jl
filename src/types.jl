#=
    Declarations of the main `types` used to define Ito diffusions
    !! IMPORTANT !!
    to define a diffusion one needs to define a struct and provide a drift
    and a volatility coefficient functions:
    `b`: denotes a drift
    `σ`: denotes a volatility coefficient
    `a`: (=:σσ') denotes a diffusivity coefficient
    for linear diffusions
    `B` and `β` jointly define the drift via `b(t,x):=Bₜx+βₜ`.
=#
using LinearAlgebra
import Base: lowercase, eltype

"""
    DiffusionProcess{T,DP,DW}

Types inheriting from `DiffusionProcess` define Ito diffusions. `T` denotes the
datatype of each coordinate, `DP` the dimension of the stochastic process,
`DW` the dimension of the Wiener process, `SS` lists the state space
restrictions.
"""
abstract type DiffusionProcess{T,DP,DW,SS} end

a(t, x, P::DiffusionProcess) = σ(t, x, P) * σ(t, x, P)'

"""
    dimension(::DiffusionProcess{T,DP,DW})

Return dimension of the stochastic process and driving Brownian motion.
"""
dimension(::DiffusionProcess{T,DP,DW}) where {T,DP,DW} = (
    process = DP,
    wiener = DW
)

"""
    eltype(::DiffusionProcess{T}) where T = T

Return the datatype that each coordinate of the stochastic process is stored in.
"""
eltype(::DiffusionProcess{T}) where T = T

"""
    state_space(::DiffusionProcess{T,DP,DW,SS})

Return the state space restrictions.
"""
state_space(::DiffusionProcess{T,DP,DW,SS}) where {T,DP,DW,SS} = SS

"""
    LinearDiffusion{T,DP,DW,SS} <: DiffusionProcess{T,DP,DW,SS}

Types inheriting from `LinearDiffusion` define a linear Ito-type diffusion, i.e.
solutions to stochastic differential equations of the form:
dXₜ = (BₜXₜ + βₜ)dt + σₜdWₜ, t∈[0,T], X₀=x₀.
"""
abstract type LinearDiffusion{T,DP,DW,SS} <: DiffusionProcess{T,DP,DW,SS} end

b(t, x, P::LinearDiffusion) = B(t, P)*x + β(t, P)
σ(t, x, P::LinearDiffusion) = σ(t, P)
a(t, x, P::LinearDiffusion) = a(t, P)
a(t, P::LinearDiffusion) = σ(t, P) * σ(t, P)'

"""
    DiffusionDomain

Types inheriting from `DiffusionStateSpace` define the types of restrictions put
on the state-space of the stochastic process.
"""
abstract type DiffusionStateSpace end
