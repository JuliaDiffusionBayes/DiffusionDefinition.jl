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
    struct Wiener{D}
    end

A struct defining the Wiener process. `D` indicates the dimension of the process
if it cannot be inferred from the DataType of the Trajectory. If the dimension
can be inferred then it takes precedence over the value of `D`.
"""
struct Wiener{D}
    Wiener{D}() where D = new{D}()
    Wiener(D::Integer) = new{D}()
    Wiener() = new{false}()
end
wiener() = Wiener()

ismutable(el) = ismutable(typeof(el))
ismutable(::Type) = Val(false)
ismutable(::Type{<:Array}) = Val(true)

"""
    Base.zero(K::Type, D, ::Val)

If `K` is a mutable type, then create `zeros` of dimension `D` and entries with
types `eltype(K)`. Otherwise, calls regular zero(`K`).
"""
Base.zero

Base.zero(K::Type, D, ::Val{true}) = zeros(eltype(K), D)
Base.zero(K::Type, D, ::Val{false}) = zero(K)

"""
    Trajectories.trajectory(tt, v::Type, D::Number)

Create a Trajectory with mutable states of dimension `D` along the time
collection `tt`.
"""
function Trajectories.trajectory(tt, v::Type, D::Number)
    trajectory(collect(tt), [zeros(eltype(v),D) for _ in tt])
end

"""
    Trajectories.trajectory(tt, v::Type)

Create a Trajectory with immutable states along the time collection `tt`.
"""
function Trajectories.trajectory(tt, v::Type)
    trajectory(collect(tt), zeros(v, length(tt)))
end
Trajectories.trajectory(tt, v::Type, D, ::Val{true}) = trajectory(tt, v, D)
Trajectories.trajectory(tt, v::Type, D, ::Val{false}) = trajectory(tt, v)

"""
    Random.rand!(
        path::Trajectory{T,Vector{K}},
        w::Wiener{D},
        y1=zero(K,D,ismutable(K))
    ) where {T,K,D}

Samples Wiener process started from `y1` and saves the data in `path`. Uses
a default random number generator.
"""
function Random.rand!(
        path::Trajectory{T,Vector{K}},
        w::Wiener{D},
        y1=zero(K,D,ismutable(K))
    ) where {T,K,D}
    rand!(Random.GLOBAL_RNG, path, w, y1)
end

"""
    Random.rand!(
        rng::Random.AbstractRNG,
        path::Trajectory{T,Vector{K}},
        ::Wiener,
        y1=zero(K)
    ) where {T,K}

Samples Wiener process with immutable states, started from `y1` and saves the
data in `path`. `rng` is used as a random number generator.
"""
function Random.rand!(
        rng::Random.AbstractRNG,
        path::Trajectory{T,Vector{K}},
        ::Wiener,
        y1=zero(K)
    ) where {T,K}
    yy, tt = path.x, path.t
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
        path::Trajectory{T,Vector{Vector{K}}},
        ::Wiener{D},
        y1=zeros(K,D)
    ) where {K,T,D}

Samples Wiener process with mutable states, started from `y1` and saves the data
in `path`. `rng` is used as a random number generator.
"""
function Random.rand!(
        rng::Random.AbstractRNG,
        path::Trajectory{T,Vector{Vector{K}}},
        ::Wiener{D},
        y1=zeros(K,D)
    ) where {K,T,D}
    yy, tt = path.x, path.t
    N = length(path)
    yy[1] = y1
    for i in 2:N
        rootdt = sqrt(tt[i] - tt[i-1])
        for j in 1:D
            yy[i][j] = yy[i-1][j] + rootdt*randn(rng, K)
        end
    end
    path
end

"""
    Base.rand(tt, w::Wiener, y1)

Samples Wiener process on `tt`, started from `y1` and returns a new object with
a sampled trajectory.
"""
Base.rand(tt, w::Wiener, y1) = rand(Random.GLOBAL_RNG, tt, w, y1)

function Base.rand(
        rng::Random.AbstractRNG,
        tt,
        w::Wiener{D},
        y1::K
    ) where {D,K}
    path = trajectory(tt, K, D, ismutable(K))
    rand!(rng, path, w, y1)
end

"""
    solve!(XX, WW, P, y1)

Compute a trajectory, started from `y1` and following the diffusion law `P`,
from the sampled Wiener process `WW`. Save the sampled path in `XX`. Return
prematurely with a `false` massage if the numerical scheme has led to the solver
violating the state-space restrictions.
"""
solve!(XX, WW, P, y1) = solve!(EulerMaruyama(), XX, WW, P, y1)

"""
    solve!(XX, WW, P, y1, buffer)

Same as `solve!(XX, WW, P, y1)`, but additionally provides a pre-allocated
buffer for performing in-place computations.
"""
solve!(XX, WW, P, y1, buffer) = solve!(EulerMaruyama(), XX, WW, P, y1, buffer)

function solve!(
        ::EulerMaruyama,
        XX::Trajectory{T,Vector{K}},
        WW::Trajectory{T,Vector{K}},
        P,
        y1::K,
    ) where {K,T}
    yy, ww, tt = XX.x, WW.x, XX.t
    N = length(XX)

    yy[1] = y1
    for i in 2:N
        dt = tt[i] - tt[i-1]
        dW = ww[i] - ww[i-1]
        yy[i] = yy[i-1] + b(tt[i-1], yy[i-1], P)*dt + σ(tt[i-1], yy[i-1], P)*dW
        bound_satisfied(P, yy[i]) || return false
    end
    true
end

function solve!(
        ::EulerMaruyama,
        XX::Trajectory{T,Vector{Vector{K}}},
        WW::Trajectory{T,Vector{Vector{K}}},
        P,
        y1::Vector{K},
        buffer=StandardEulerBuffer{K}(dimensions(P)...)
    ) where {K,T}
    yy, ww, tt = XX.x, WW.x, XX.t
    N = length(XX)
    dim_proc = dimension(P).process

    yy[1] = y1
    dt = tt[2] - tt[1]
    y = yy[1]

    for i in 2:N
        dt = tt[i] - tt[i-1]
        y = yy[i-1]
        b!(buffer, tt[i-1], y, P)
        σ!(buffer, tt[i-1], y, P)
        for j in 1:dim_proc
            buffer.dW[j] = ww[i][j] - ww[i-1][j]
        end
        mul!(yy[i], buffer.σ, buffer.dW)
        for j in 1:dim_proc
            yy[i][j] += y[j] + buffer.b[j]*dt
        end
        bound_satisfied(P, yy[i]) || return false
    end
    true
end
