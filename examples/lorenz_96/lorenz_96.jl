@diffusion_process Lorenz96 begin
    :dimensions
    process --> DiffusionDefinition._TEMP_DIMENSION_PROCESS
    wiener --> DiffusionDefinition._TEMP_DIMENSION_WIENER

    :parameters
    (θ, σ) --> Float64

    :additional
    constdiff --> true
    diagonaldiff --> true
end

#-------------------------------------------------------------------------------
#                Some helper functions for static indexing

function si1(::DiffusionDefinition.DiffusionProcess{T,DP}) where {T,DP}
    SVector{DP,Int64}(mod1.(2:DP+1, DP))
end

function si2(::DiffusionDefinition.DiffusionProcess{T,DP}) where {T,DP}
    SVector{DP,Int64}(mod1.(-1:DP-2, DP))
end

function si3(::DiffusionDefinition.DiffusionProcess{T,DP}) where {T,DP}
    SVector{DP,Int64}(mod1.(0:DP-1, DP))
end

function si4(::DiffusionDefinition.DiffusionProcess{T,DP}) where {T,DP}
    SVector{DP,Int64}(1:DP)
end

function sid_mat(::DiffusionDefinition.DiffusionProcess{T,DP,DW}) where {T,DP,DW}
    SDiagonal{DW,T}(I)
end

#-------------------------------------------------------------------------------
#                Some helper functions for dynamic indexing

function i1(::DiffusionDefinition.DiffusionProcess{T,DP}) where {T,DP}
    mod1.(2:DP+1, DP)
end

function i2(::DiffusionDefinition.DiffusionProcess{T,DP}) where {T,DP}
    mod1.(-1:DP-2, DP)
end

function i3(::DiffusionDefinition.DiffusionProcess{T,DP}) where {T,DP}
    mod1.(0:DP-1, DP)
end

function i4(::DiffusionDefinition.DiffusionProcess{T,DP}) where {T,DP}
    1:DP
end

function id_mat(::DiffusionDefinition.DiffusionProcess{T,DP,DW}) where {T,DP,DW}
    Diagonal{T}(I, DW)
end

#-------------------------------------------------------------------------------
#                 Drift and diffusion coefficients

function b(t, x, P::Lorenz96)
    (
        (view(x, si1(P)) - view(x, si2(P))) .* view(x, si3(P))
        - view(x, si4(P)) .+ P.θ
    )
end

σ(t, x, P::Lorenz96) = P.σ * sid_mat(P)

function b!(buffer, t, x, P::Lorenz96)
    buffer.b .= (
        (view(x, i1(P)) .- view(x, i2(P))) .* view(x, i3(P))
        .- view(x, i4(P))
        .+ P.θ
    )
end

σ!(buffer, t, x, P::Lorenz96) = (buffer.σ.diag .= P.σ)
