@diffusion_process Lorenz{T} begin
    :dimensions
    process --> 3
    wiener --> 3

    :parameters
    _ --> (3, T)
    σ --> T

    :additional
    constdiff --> true
    sparsediff --> true
end

function b(t, x, P::Lorenz)
    @SVector [
        P.p1*(x[2]-x[1]),
        P.p2*x[1] - x[2] - x[1]*x[3],
        x[1]*x[2] - P.p3*x[3]
    ]
end

σ(t, x, P::Lorenz) = SDiagonal(P.σ, P.σ, P.σ)
