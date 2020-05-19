@diffusion_process JansenRit{T} begin
    :dimensions
    process --> 6
    wiener --> 1
    :parameters
    (A, a, B, b, C, α1, α2, νmax, v, r, μy, σy) --> T
    :additional
    constdiff --> true
end
# θ
#param_names(P::JansenRit) = (:A, :a, :B, :b, :C, :α1, :α2, :νmax, :v, :r, :μy, :σy)
#θ =[3.25, 100.0, 22.0, 50.0, 135.0, 0.8, 0.25, 5.0, 6.0, 0.56, 220.0, 2000.0]
#P = JansenRit(θ...)

sigm(x, P::JansenRit) = P.νmax / (1 + exp(P.r*(P.v - x)))
μy(t, P::JansenRit) = P.μy #constant
C1(P::JansenRit) = P.C
C2(P::JansenRit) = P.α1*P.C
C3(P::JansenRit) = P.α2*P.C
C4(P::JansenRit) = P.α2*P.C

function b(t, x, P::JansenRit)
    @SVector [
        x[4],
        x[5],
        x[6],
        P.A*P.a*(sigm(x[2] - x[3], P)) - 2P.a*x[4] - P.a*P.a*x[1],
        P.A*P.a*(μy(t, P) + C2(P)*sigm(C1(P)*x[1], P)) - 2P.a*x[5] - P.a*P.a*x[2],
        P.B*P.b*(C4(P)*sigm(C3(P)*x[1], P)) - 2P.b*x[6] - P.b*P.b*x[3]]
end

function σ(t, x, P::JansenRit)
    @SMatrix [
        0.0;
        0.0;
        0.0;
        0.0;
        P.σy;
        0.0
    ]
end
