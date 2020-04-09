#==============================================================================#
#
#    Declarations of the main `types` used to define Ito diffusions
#    !! IMPORTANT !!
#    to define a diffusion one needs to define a struct and provide a drift
#    and a volatility coefficient functions:
#    `b`: denotes a drift
#    `σ`: denotes a volatility coefficient
#    `a`: (=:σσ') denotes a diffusivity coefficient
#    for linear diffusions
#    `B` and `β` jointly define the drift via `b(t,x):=Bₜx+βₜ`.
#
#==============================================================================#

"""
    DiffusionProcess{T,DP,DW,SS,EI}

Types inheriting from `DiffusionProcess` define Ito diffusions. `T` denotes the
datatype of each coordinate, `DP` the dimension of the stochastic process,
`DW` the dimension of the Wiener process, `SS` lists the state space
restrictions
"""
abstract type DiffusionProcess{T,DP,DW,SS} end

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


const IndexedTime = Tuple{T,S} where {T<:Number,S<:Integer}

#------------------------------------------------------------------------------#
#       Names for functions defining drifts and volatility coefficients
#------------------------------------------------------------------------------#

"""
    b

Compute the drift of a diffusion out-of-place (should use StaticArrays).
"""
function b end

"""
    B

Compute matrix `B` of a linear diffusion out-of-place (should use StaticArrays).
"""
function B end

"""
    β

Compute vector `β` of a linear diffusion out-of-place (should use StaticArrays).
"""
function β end

"""
    σ

Compute the volatility coefficient of a diffusion out-of-place (should use
StaticArrays).
"""
function σ end

"""
    a

Compute the diffusion function σσ' of a diffusion out-of-place (should use
StaticArrays).
"""
function a end

"""
    b!

Compute the drift of a diffusion in-place (uses buffers to store temporary
results).
"""
function b! end

"""
    B!

Compute matrix `B` of a linear diffusion in-place (uses buffers to store
temporary results).
"""
function B! end

"""
    β!

Compute vector `β!` of a linear diffusion in-place (uses buffers to store
temporary results).
"""
function β! end

"""
    σ!

Compute the volatility coefficient of a diffusion in-place (uses buffers to
store temporary results).
"""
function σ! end

"""
    a!

Compute the diffusion function σσ' of a diffusion in-place (uses buffers to
store temporary results).
"""
function a! end

"""
    constdiff

Returns true if the diffusion coefficient does not depend on the state of the
process or time.
"""
function constdiff end
