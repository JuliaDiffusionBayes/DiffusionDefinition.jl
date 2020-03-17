@diffusion_process JansenRit begin
    :dimensions
    process --> 6
    wiener --> 3

    :parameters
    (A, a, B, b, C, νmax, v, r, μx, μy, μz, σx, σy, σz) --> Float64

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

    :additional
    constdiff --> true
end

sigm(x, P::JansenRit) = P.νmax / (1 + exp(P.r*(P.v0 - x)))
μx(t, P::JansenRit) = P.μx
μy(t, P::JansenRit) = P.μy
μz(t, P::JansenRit) = P.μz

function b(t, x, P::JansenRit)
    ℝ{6}(
        x[4],
        x[5],
        x[6],
        P.δ*x[1]*x[2] - P.γ*x[2],
        P.A*P.a*(μx(t, P) + sigm(x[2] - x[3], P)) - 2P.a*x[4] - P.a*P.a*x[1],
        P.A*P.a*(μy(t, P) + 0.8P.C*sigm(P.C*x[1], P)) - 2P.a*x[5] - P.a*P.a*x[2],
        P.B*P.b*(μz(t, P) + 0.25P.C*sigm(0.25P.C*x[1], P)) - 2P.b*x[6] - P.b*P.b*x[3])
    )
end

function σ(t, x, P::JansenRit)
    @SMatrix [
        0.0  0.0  0.0;
        0.0  0.0  0.0;
        0.0 0.0  0.0;
        P.σx  0.0  0.0;
        0.0  P.σy 0.0;
        0.0  0.0  P.σz
    ]
end
