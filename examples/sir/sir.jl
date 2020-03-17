@diffusion_process SIR begin
    :parameters
    (α, β, σ1, σ2) --> Float64

    :conjugate
    phi(t, u) --> (
        (0.0, 0.0),
        ((k - u[1] - u[2])*u[1], 0.0),
        (-u[1], u[1]),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
    )
    nonhypo(x) --> x
    num_non_hypo --> 2

    :additional
    statespace -->BoundedStateSpace(((1,2), (0.0,0.0)), ((1,2), (1.0,1.0)))
end

function b(t, x, P::SIR)
    ℝ{2}(
        P.α*(1 - u[1] - u[2])*u[1] - P.β*u[1],
        P.β*u[1]
    )
end

function σ(t, x, P::SIR)
    @SMatrix Float64[
    -P.σ1*sq((1 - u[1] - u[2])*u[1])  -P.σ2*sq(u[1]);
    0.0   P.σ2*sq.(u[1])
]
