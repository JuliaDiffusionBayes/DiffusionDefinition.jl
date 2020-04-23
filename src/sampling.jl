#==============================================================================#
#
#               Routines for sampling trajectories from the
#               diffusion laws defined by this package.
#
#==============================================================================#
"""
    abstract type AbstractSDESolver end

Supertype of flags indicating ODE solvers
"""
abstract type AbstractSDESolver end

"""
    struct EulerMaruyama <: AbstractSDESolver end

Flag for indicating use of the Euler-Maruyama scheme for sampling diffusions.
IMPORTANT: this is the only diffusion path sampler implemented in this package.
There are no plans for implementing any other SDE solver in the forseeable
future.
"""
struct EulerMaruyama <: AbstractSDESolver end

"""
    struct Wiener{D,Tdevice}
    end

A struct defining the Wiener process. `D` indicates the dimension of the process
if it cannot be inferred from the DataType of the Trajectory. If the dimension
can be inferred then it takes precedence over the value of `D`.
"""
struct Wiener{D,Tdevice}
    Wiener{D,Tdevice}() where {D,Tdevice} = new{D,Tdevice}()
    Wiener{D}() where D = new{D,:unspec}()
    Wiener(D::Integer=false, device::Symbol=:unspec) = new{D,device}()
    function Wiener(
            ::DiffusionProcess{T,TP,TW},
            device::Symbol=:unspec
        ) where {T,TP,TW}
        new{TW,device}()
    end
end
wiener(args...) = Wiener(args...)

ismutable(el) = ismutable(typeof(el))
ismutable(::Type) = Val(true)
ismutable(::Type{K}) where K<:Union{Number,SArray} = Val(false)
ismutable(::Trajectory{T,Vector{K}}) where {T,K} = ismutable(K)

Base.zero(w::Wiener{D}) where D = zero(SVector{D,Float64})

#===============================================================================
                        Sampling Wiener processes
===============================================================================#
"""
    Base.rand(w::Wiener, tt, y1)

Samples Wiener process on `tt`, started from `y1` and returns a new object with
a sampled trajectory.
"""
Base.rand(w::Wiener, tt, y1=zero(w)) = rand(Random.GLOBAL_RNG, w, tt, y1)

function Base.rand(
        rng::Random.AbstractRNG,
        w::Wiener{D},
        tt,
        y1::K=zero(w),
    ) where {D,K}
    path = trajectory(tt, K, D, ismutable(K))
    rand!(rng, w, path, y1)
end

"""
    Random.rand!(
        rng::Random.AbstractRNG,
        ::Wiener,
        path::Trajectory{T,Vector{K}},
        y1=zero(K)
    ) where {T,K}

Samples Wiener process with immutable states, started from `y1` and saves the
data in `path`. `rng` is used as a random number generator.
"""
function Random.rand!(
        rng::Random.AbstractRNG,
        ::Wiener,
        path::Trajectory{T,<:Vector{K}},
        y1=zero(K)
    ) where {T,K}
    tt, yy = path.t, path.x

    N = length(path)
    yy[1] = y1
    for i in 2:N
        rootdt = sqrt(tt[i] - tt[i-1])
        yy[i] = yy[i-1] + rootdt*randn(rng, K)
    end
    path
end

"""
    Random.rand!(
        rng::Random.AbstractRNG,
        ::Wiener{D},
        path::Trajectory{T,Vector{Vector{K}}},
        y1=zeros(K,D)
    ) where {K,T,D}

Samples Wiener process with mutable states, started from `y1` and saves the data
in `path`. `rng` is used as a random number generator.
"""
function Random.rand!(
        rng::Random.AbstractRNG,
        ::Wiener{D},
        path::Trajectory{T,<:Vector{<:Array{K}}},
        y1=zeros(K,D)
    ) where {K,T,D}
    yy, tt = path.x, path.t
    N = length(path)
    yy[1] = y1
    for i in 2:N
        rootdt = sqrt(tt[i] - tt[i-1])
        for j in eachindex(yy[i])
            yy[i][j] = yy[i-1][j] + rootdt*randn(rng, K)
        end
    end
    path
end

"""
    Random.rand!(
        w::Wiener{D},
        path::Trajectory{T,Vector{K}},
        y1=zero(K,D,ismutable(K))
    ) where {T,K,D}

Samples Wiener process started from `y1` and saves the data in `path`. Uses
a default random number generator.
"""
function Random.rand!(
        w::Wiener{D},
        path::Trajectory{T,Vector{K}},
        y1=zero(K, D, ismutable(K))
    ) where {T,K,D}
    rand!(Random.GLOBAL_RNG, w, path, y1)
end

#===============================================================================
                    Euler-Maruyama scheme for solving SDEs
===============================================================================#

# needed for ForwardDiff
value(x) = x
value!(y, x) = (y.= x)

"""
    solve!(XX, WW, P, y1)

Compute a trajectory, started from `y1` and following the diffusion law `P`,
from the sampled Wiener process `WW`. Save the sampled path in `XX`. Return
prematurely with a `false` massage if the numerical scheme has led to the solver
violating the state-space restrictions.
"""
solve!(XX, WW, P, y1) = solve!(EulerMaruyama(), ismutable(XX), XX, WW, P, y1)

"""
    solve!(XX, WW, P, y1, buffer)

Same as `solve!(XX, WW, P, y1)`, but additionally provides a pre-allocated
buffer for performing in-place computations.
"""
function solve!(XX, WW, P, y1, buffer)
    solve!(EulerMaruyama(), ismutable(XX), XX, WW, P, y1, buffer)
end

function solve!(::EulerMaruyama, ::Val{false}, XX, WW, P, y1::Ky1) where Ky1
    yy, ww, tt = XX.x, WW.x, XX.t
    N = length(XX)

    yy[1] = value(y1)
    y = y1 # more appropriate name
    for i in 2:N
        dt = tt[i] - tt[i-1]
        dW = ww[i] - ww[i-1]
        y = y + _b((tt[i-1], i-1), y, P)*dt + _σ((tt[i-1], i-1), y, P)*dW
        yy[i] = value(y) # strip duals
        bound_satisfied(P, yy[i]) || return false, nothing
    end
    true, y
end

function solve!(
        ::EulerMaruyama, ::Val{true}, XX, WW, P, y1::K,
        buffer=StandardEulerBuffer{K}(P)
    ) where K
    yy, ww, tt = XX.x, WW.x, XX.t
    N = length(XX)

    value!(yy[1], y1)
    y::K = y1 # name alias

    for i in 2:N
        dt = tt[i] - tt[i-1]

        _b!(buffer, (tt[i-1], i-1), y, P)
        _σ!(buffer, (tt[i-1], i-1), y, P)

        buffer.dW .= ww[i] .- ww[i-1]
        mul!(buffer.y, buffer.σ, buffer.dW)
        y .= y .+ buffer.y .+ buffer.b .* dt

        value!(yy[i], y) # strip duals
        bound_satisfied(P, yy[i]) || return false, nothing
    end
    true, y
end


#=
function solve!(
        ::EulerMaruyama, ::Val{true}, XX, WW, P, y1::K,
        buffer=StandardEulerBuffer{K}(P)
    ) where K
    yy, ww, tt = XX.x, WW.x, XX.t
    N = length(XX)

    value!(yy[1], y1)
    y = y1 # name alias
    y_temp = get_y(buffer)
    b_temp = get_b(buffer)
    σ_temp = get_σ(buffer)
    dW = get_dW(buffer)

    for i in 2:N
        dt = tt[i] - tt[i-1]

        _b!(buffer, (tt[i-1], i-1), y, P)
        _σ!(buffer, (tt[i-1], i-1), y, P)

        dW .= ww[i] .- ww[i-1]
        mul!(y_temp, σ_temp, dW)
        y .+= y_temp .+ b_temp .* dt

        value!(yy[i], y) # strip duals
        bound_satisfied(P, yy[i]) || return false, nothing
    end
    true, y
end
=#

#===============================================================================
        Sampling Diffusion processes using the Euler Maruyama scheme
===============================================================================#

function Base.rand(P::DiffusionProcess, tt, y1=zero(P))
    rand(Random.GLOBAL_RNG, P, tt, y1, ismutable(y1))
end

function Base.rand(
        rng::Random.AbstractRNG,
        P::DiffusionProcess{T,DP,DW},
        tt,
        y1::K,
        v::Val{true}
    ) where {T,DP,DW,K}
    WW = rand(rng, Wiener{DW}(), tt, zeros(eltype(y1), DW))
    XX = trajectory(tt, K, DP, v)
    success, end_pt = false, nothing
    buffer = StandardEulerBuffer{K}(P)
    while !success
        success, end_pt = solve!(XX, WW, P, y1, buffer)
    end
    XX, end_pt
end

StaticArrays.similar_type(::Type{K}, s::Size) where K <: Number = SVector{s[1],K}

function Base.rand(
        rng::Random.AbstractRNG,
        P::DiffusionProcess{T,DP,DW},
        tt,
        y1::K,
        v::Val{false}
    ) where {T,DP,DW,K}
    WW = rand(rng, Wiener(), tt, zero(similar_type(K, Size(DW))))
    XX = trajectory(tt, K, DP, v)
    success, end_pt = false, nothing
    while !success
        success, end_pt = solve!(XX, WW, P, y1)
    end
    XX, end_pt
end
