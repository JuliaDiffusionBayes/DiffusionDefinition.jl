#=
    Declarations of the main `types` used to define Ito diffusions
=#

import Base: lowercase, eltype

"""
    DiffusionProcess{T,DP,DW}

Types inheriting from `DiffusionProcess` define Ito diffusions. `T` denotes the
datatype of each coordinate, `DP` the dimension of the stochastic process,
`DW` the dimension of the Wiener process, `SS` lists the state space
restrictions.
"""
abstract type DiffusionProcess{T,DP,DW,SS} end

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

"""
    DiffusionDomain

Types inheriting from `DiffusionStateSpace` define the types of restrictions put
on the state-space of the stochastic process.
"""
abstract type DiffusionStateSpace end
