#NOTE Engines for random numbers need to be figured out, currently falls on
# defaults
using CuArrays
using CUDAnative
using Random

const _GLOBAL_NUM_THREADS = 256

Base.zero(w::Wiener{D,:GPU}) where D = CuArrays.fill(0.0f0, D)

struct GPUTrajectory{S,T} <: AbstractPairedArray
    t::S
    x::T
    dt::S
    sqdt::S
end

function gpu_trajectory(t::S, x::T) where {S,T}
    GPUTrajectory{S,T}(t, x, t[2:end] .- t[1:end-1], sqrt.(t[2:end] .- t[1:end-1]))
end

#===============================================================================
                        Sampling Wiener processes
===============================================================================#

function Base.rand(w::Wiener{D,:GPU}, tt, y1::CuArray{K}=zero(w)) where {D,K}
    N = length(y1)
    path = gpu_trajectory(tt, zeros(CuArray{K}, (N, length(tt))); sqdt=true)
    rand!(w, path.tr, y1)
end

function Random.rand!(
        ::Wiener{D},
        path::GPUTrajectory,
        y1=zero(path.x[:,1])
    ) where D
    numblocks = ceil(Int, D/_GLOBAL_NUM_THREADS)

    yy = path.x
    randn!(yy)
    yy[:,1] .= y1
    @cuda threads=_GLOBAL_NUM_THREADS blocks=numblocks solve_Wiener_GPU!(yy, path.sqdt)
    path
end

function solve_Wiener_GPU!(y, sqdt)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = 2:size(y,2)
        for j = index:stride:size(y,1)
            @inbounds y[j,i] = y[j,i-1] + sqdt[i-1] * y[j,i]
        end
        sync_threads()
    end
    return
end


#===============================================================================
                    Euler-Maruyama scheme for solving SDEs
===============================================================================#

function solve!(
        ::EulerMaruyama, XX::GPUTrajectory, WW::GPUTrajectory, P, y1::K,
        buffer=StandardEulerBuffer{K}(P)
    ) where K
    yy, ww, tt = XX.x, WW.x, XX.t
    N = length(XX)

    value!(yy[:,1], y1)
    y = y1 # name alias

    for i in 2:N
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

function solve_Euler_GPU!(y, w, sqdt)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = 2:size(y,2)
        _b!(buffer, (tt[i-1], i-1), y, P)


        for j = index:stride:size(y,1)
            @inbounds y[j,i] = y[j,i-1] + sqdt[i-1] * y[j,i]
        end
        sync_threads()
    end
end



#===============================================================================
        Sampling Diffusion processes using the Euler Maruyama scheme
===============================================================================#

function Base.rand(
        P::DiffusionDefinition.DiffusionProcess{T,DP,DW},
        tt::CuArray,
        y1::K=zero(P)
    ) where {T,DP,DW,K<:CuArray}
    WW = rand(Wiener{DW}(), tt, zeros(K, DW))
    XX = trajectory(tt, CuArray[zeros(K, DP) for _ in tt])
    success = false
    buffer = DiffusionDefinition.StandardEulerBuffer{K}(P)
    while !success
        success = DiffusionDefinition.solve!(XX, WW, P, y1, buffer)
    end
    XX
end
