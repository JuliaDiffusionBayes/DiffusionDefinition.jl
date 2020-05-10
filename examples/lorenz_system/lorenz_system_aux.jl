@diffusion_process LorenzAux{K,R} begin
    :dimensions
    process --> 3
    wiener --> 3

    :parameters
    _ --> (3, K)
    σ --> K

    :auxiliary_info
    t0 --> Float64
    T --> Float64
    vT --> R

    :additional
    constdiff --> true
    linear --> true
end

DiffusionDefinition.B(t, P::LorenzAux) = @SMatrix [
    -P.p1         P.p1        0.0;
     P.p2-P.vT[3]  -1.0    -P.vT[1];
     P.vT[2]      P.vT[1]    -P.p3
]

DiffusionDefinition.β(t, P::LorenzAux) = @SVector[
    0.0,
    P.vT[1]*P.vT[3],
    -P.vT[1]*P.vT[3]
]

DiffusionDefinition.σ(t, P::LorenzAux) = SDiagonal(P.σ, P.σ, P.σ)
