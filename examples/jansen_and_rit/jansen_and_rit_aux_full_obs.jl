@diffusion_process JansenRitAuxFull{K, R} begin
    :dimensions
    process --> 6
    wiener --> 1
    :parameters
    (A, a, B, b, C, α1, α2, νmax, v, r, μy, σy) --> K
    :auxiliary_info
    t0 --> Float64
    T --> Float64
    vT --> R
    :additional
    constdiff --> true
    linear --> true
end

sigm(x, P::JansenRitAuxFull) = P.νmax / (1 + exp(P.r*(P.v - x)))
C1(P::JansenRitAuxFull) = P.C
C2(P::JansenRitAuxFull) = P.α1*P.C
C3(P::JansenRitAuxFull) = P.α2*P.C
C4(P::JansenRitAuxFull) = P.α2*P.C
sigm1d(x, P::JansenRitAuxFull) = P.νmax*P.r*exp(P.r*(P.v-x))/(1 + exp(P.r*(P.v-x)))^2

function B(t, P::JansenRitAuxFull)
    @SMatrix [0.0 0.0 0.0 1.0 0.0 0.0;
                0.0 0.0 0.0 0.0 1.0 0.0;
                0.0 0.0 0.0 0.0 0.0 1.0;
                -P.a^2 P.A*P.a*sigm1d(P.vT[2]-P.vT[3], P) -P.A*P.a*sigm1d(P.vT[2]-P.vT[3], P) -2.0*P.a 0.0 0.0;
                P.A*P.a*C2(P)*sigm1d(C1(P)*P.vT[1], P)*C1(P) -P.a^2 0.0 0.0 -2.0*P.a 0.0;
                P.B*P.b*C4(P)*sigm1d(C3(P)*P.vT[1], P)*C3(P) 0.0 -P.b^2 0.0 0.0 -2.0*P.b]
end



# β(t, P::JansenRitAuxFull) = @SVector [0.0, 0.0, 0.0, -P.A*P.a*sigm1d(P.vT[1], P)*P.vT[1],
                                # P.A*P.a*(P.μy + C2(P)*sigm(C1(P)*P.x1, P)),  P.B*P.b*C4(P)*sigm(C3(P)*P.x1, P)]

β(t, P::JansenRitAuxFull) = @SVector [0.0, 0.0, 0.0,
                P.A*P.a*(sigm(P.vT[2] - P.vT[3], P) - sigm1d(P.vT[2] - P.vT[3], P)*(P.vT[2] - P.vT[3])),
                P.A*P.a*(P.μy + C2(P)*(sigm(C1(P)*P.vT[1], P) - sigm1d(C1(P)*P.vT[1], P)*C1(P)*P.vT[1])),
                P.B*P.b*(C4(P)*(sigm(C3(P)*P.vT[1], P) - sigm1d(C3(P)*P.vT[1], P)*C3(P)*P.vT[1])) ]



function σ(t, P::JansenRitAuxFull)
    @SVector [0.0, 0.0, 0.0, 0.0, P.σy, 0.0]
end
