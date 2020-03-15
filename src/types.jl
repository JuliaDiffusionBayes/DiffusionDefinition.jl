#=
    Declarations of the main `types` used to define Ito diffusions
=#

import Base: lowercase, eltype

"""
    DiffusionProcess{T,DP,DW}

Types inheriting from `DiffusionProcess` define Ito diffusions. `T` denotes the
datatype of each coordinate, `DP` the dimension of the stochastic process,
`DW` the dimension of the Wiener process.
"""
abstract type DiffusionProcess{T,DP,DW} end

"""
    dimension(d::DiffusionProcess{T,DP,DW})

Return dimension of the stochastic process and driving Brownian motion.
"""
dimension(d::DiffusionProcess{T,DP,DW}) where {T,DP,DW} = (
    process = DP,
    wiener = DW
)

"""
    eltype(d::DiffusionProcess{T}) where T = T

Return the datatype that each coordinate of the stochastic process is stored in.
"""
eltype(d::DiffusionProcess{T}) where T = T

"""
    LinearDiffusion{T,DP,DW} <: DiffusionProcess{T,DP,DW}

Types inheriting from `LinearDiffusion` define a linear Ito-type diffusion, i.e.
solutions to stochastic differential equations of the form:
dXₜ = (BₜXₜ + βₜ)dt + σₜdWₜ, t∈[0,T], X₀=x₀.
"""
abstract type LinearDiffusion{T,DP,DW} <: DiffusionProcess{T,DP,DW} end

"""
    DiffusionDomain

Types inheriting from `DiffusionStateSpace` define the types of restrictions put
on the state-space of the stochastic process.
"""
abstract type DiffusionStateSpace end
