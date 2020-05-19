@diffusion_process JansenRitAux{K, R} begin
    :dimensions
    process --> 6
    wiener --> 1
    :parameters
    (a, b, σy) --> K
    :auxiliary_info
    t0 --> Float64
    T --> Float64
    vT --> R

    :additional
    constdiff --> true
    linear --> true
end
# θ
#param_names(P::JansenRitAux) = (:a, :b,, :σy)
#θ =[100.0, 50.0, 2000.0]
#P = JansenRitAux(θ...)


function B(t, P::JansenRitAux)
    @SMatrix [  0.0 0.0 0.0 1.0 0.0 0.0;
                0.0 0.0 0.0 0.0 1.0 0.0;
                0.0 0.0 0.0 0.0 0.0 1.0;
                -2.0*P.a 0.0 0.0 -P.a^2 0.0 0.0;
                0.0 -2.0*P.a 0.0 0.0 -P.a^2 0.0;
                0.0 0.0 -2.0*P.b 0.0 0.0 -P.b^2]
end


β(t, P::JansenRitAux) = @SVector [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


function σ(t, P::JansenRitAux)
    @SVector [0.0, 0.0, 0.0, 0.0, P.σy, 0.0]
end
