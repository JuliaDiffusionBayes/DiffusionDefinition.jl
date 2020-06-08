@diffusion_process LotkaVolterra{K} begin
    :dimensions
    process --> 2
    wiener --> 2

    :parameters
    (α, β, γ, δ, σ1, σ2) --> K

    #=
    :conjugate
    phi(t, x) --> (
        (0.0, 0.0),
        (x[1], 0.0),
        (-x[1]*x[2], 0.0),
        (0.0, x[1]*x[2]),
        (0.0, -x[2]),
        (0.0, 0.0),
        (0.0, 0.0)
    )
    nonhypo(x) --> x
    num_non_hypo --> 2
    =#

    :additional
    constdiff --> true
    statespace --> LowerBoundedStateSpace((1,2),(0.0, 0.0))
end

function DiffusionDefinition.b(t, x, P::LotkaVolterra)
    @SVector [
        P.α*x[1] - P.β*x[1]*x[2],
        P.δ*x[1]*x[2] - P.γ*x[2],
    ]
end

DiffusionDefinition.σ(t, x, P::LotkaVolterra) = SDiagonal(P.σ1, P.σ2)
