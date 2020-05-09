@diffusion_process SIR{T} begin
    :dimensions
    process --> 2
    wiener --> 2

    :parameters
    (α, β, σ1, σ2) --> T

    #=
    :conjugate
    phi(t, u) --> (
        (0.0, 0.0),
        ((1.0 - u[1] - u[2])*u[1], 0.0),
        (-u[1], u[1]),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
    )
    nonhypo(x) --> x
    num_non_hypo --> 2
    =#

    :additional
    statespace -->BoundedStateSpace(((1,2), (0.0,0.0)), ((1,2), (1.0,1.0)))
end

function DiffusionDefinition.b(t, x, P::SIR)
    @SVector[
        P.α*(1.0 - x[1] - x[2])*x[1] - P.β*x[1],
        P.β*x[1]
    ]
end

function DiffusionDefinition.σ(t, x, P::SIR)
    @SMatrix [
        -P.σ1*sqrt((1.0 - x[1] - x[2])*x[1])  -P.σ2*sqrt(x[1]);
        0.0   P.σ2*sqrt(x[1])
    ]
end
