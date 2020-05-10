@diffusion_process JansenRit{T} begin
    :dimensions
    process --> 6
    wiener --> 3

    :parameters
    (A, a, B, b, C, νmax, ν0, r, μx, μy, μz, σx, σy, σz) --> T

    #=
    :conjugate
    phi(t, x) --> (
        (P.A*P.a*0.8P.C*sigm(P.C*x[1], P) - 2P.a*x[5] - P.a*P.a*x[2],),
        (P.A*P.a,),
        (0,),
        (0,),
        (0,),
        (0,),
        (0,),
        (0,),
        (0,),
        (0,),
        (0,),
        (0,),
        (0,),
    )
    =#

    :additional
    constdiff --> true
end

sigm(x, P::JansenRit) = P.νmax / (1 + exp(P.r*(P.ν0 - x)))
μx(t, P::JansenRit) = P.μx
μy(t, P::JansenRit) = P.μy
μz(t, P::JansenRit) = P.μz

DiffusionDefinition.b(t, x, P::JansenRit) = @SVector [
    x[4],
    x[5],
    x[6],
    P.A*P.a*(μx(t, P) + sigm(x[2] - x[3], P)) - 2*P.a*x[4] - P.a^2*x[1],
    P.A*P.a*(μy(t, P) + 0.8*P.C*sigm(P.C*x[1], P)) - 2*P.a*x[5] - P.a^2*x[2],
    P.B*P.b*(μz(t, P) + 0.25*P.C*sigm(0.25*P.C*x[1], P)) - 2*P.b*x[6] - P.b^2*x[3]
]


DiffusionDefinition.σ(t, x, P::JansenRit) = @SMatrix [
    0.0     0.0     0.0;
    0.0     0.0     0.0;
    0.0     0.0     0.0;
    P.σx    0.0     0.0;
    0.0     P.σy    0.0;
    0.0     0.0     P.σz
]
