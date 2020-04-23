#NOTE Engines for random numbers need to be figured out, currently falls on
# defaults

Base.zero(w::Wiener{D,:GPU}) where D = CuArrays.fill(0.0f0, D)

#===============================================================================
                        Sampling Wiener processes
===============================================================================#

function Base.rand(w::Wiener{D,:GPU}, tt::CuArray, y1::CuArray{K}=zero(w)) where {D,K}
    N = length(y1)
    path = trajectory(tt, CuArray[zeros(CuArray{K}, N) for _ in tt])
    rand!(w, path, y1)
end

function Random.rand!(
        ::Wiener,
        path::Trajectory{<:CuArray,<:CuArray{K}},
        y1=zeros(K, length(path.x[1]))
    ) where K
    yy, tt = path.x, path.t
    N = length(path)
    yy[1] = y1
    for i in 2:N
        rootdt = sqrt(tt[i] - tt[i-1])
        randn!(yy[i])
        yy[i] .*= rootdt
        yy[i] .+= yy[i-1]
    end
    path
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
