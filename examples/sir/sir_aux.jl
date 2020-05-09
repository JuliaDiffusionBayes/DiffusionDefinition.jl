@diffusion_process SIRAux{K,R} begin
    :dimensions
    process --> 2
    wiener --> 2

    :parameters
    (α, β, σ1, σ2) --> K

    :auxiliary_info
    t0 --> Float64
    T --> Float64
    vT --> R

    :additional
    linear --> true
    constdiff --> true
    statespace -->BoundedStateSpace(((1,2), (0.0,0.0)), ((1,2), (1.0,1.0)))
end

function B(t, P::SIRAux)
    @SMatrix[
        (P.α*(1.0 - P.vT[1] - P.vT[2]) - P.β)  0.0;
        P.β   0.0
    ]
end

β(t, P::SIRAux) = @SVector [0.0, 0.0]


function σ(t, P::SIRAux)
    @SMatrix [
        (-P.σ1*(1 - P.vT[1] - P.vT[2])*P.vT[1])  -P.σ2*P.vT[1];
        0.0   P.σ2*P.vT[1]
    ]
end
