@diffusion_process JansenRitAux2{K, R} begin
    :dimensions
    process --> 6
    wiener --> 1
    :parameters
    (A, a, B, b, νmax, v, r, μy, σy)--> K
    :auxiliary_info
    t0 --> Float64
    T --> Float64
    vT --> R
    :additional
    constdiff --> true
    linear --> true
end

function B(t, P::JansenRitAux2)
    @SMatrix [  0.0 0.0 0.0 1.0 0.0 0.0;
                0.0 0.0 0.0 0.0 1.0 0.0;
                0.0 0.0 0.0 0.0 0.0 1.0;
                -P.a^2 P.A*P.a*sigm1d(P.vT[1], P) -P.A*P.a*sigm1d(P.vT[1], P) -2.0*P.a 0.0 0.0;
                0.0 -P.a^2 0.0 0.0 -2.0*P.a 0.0;
                0.0 0.0 -P.b^2 0.0 0.0 -2.0*P.b]
end

sigm(x, P::JansenRitAux2) = P.νmax / (1 + exp(P.r*(P.v - x)))
sigm1d(x, P::JansenRitAux2) = P.νmax*P.r*exp(P.r*(P.v-x))/(1 + exp(P.r*(P.v-x)))^2
β(t, P::JansenRitAux2) = @SVector [0.0, 0.0, 0.0, P.A*P.a*(sigm(P.vT[1], P) - sigm1d(P.vT[1], P)*P.vT[1]),
                                                                P.A*P.a*P.μy,  0.0]



function σ(t, P::JansenRitAux2)
    @SVector [0.0, 0.0, 0.0, 0.0, P.σy, 0.0]
end
