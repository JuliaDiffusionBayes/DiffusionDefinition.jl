@diffusion_process SIRAux begin
    :parameters
    (α, β, σ1, σ2, t, T) --> Float64
    (u, v) --> ℝ{2}

    :additional
    linear --> true
    constdiff --> true
    statespace -->BoundedStateSpace(((1,2), (0.0,0.0)), ((1,2), (1.0,1.0)))
end

function B(t, P::SIRAux)
    @SMatrix[
        (P.α*(1 - P.v[1] - P.v[2]) - P.β)  0.0;
        P.β   0.0
    ]
end

β(t, P::SIRAux) = ℝ{2}(0.0, 0.0)


function σ(t, P::SIRAux)
    @SMatrix Float64[
        (-P.σ1*(1 - v[1] - v[2])*v[1])  -P.σ2*v[1];
        0.0   P.σ2*v[1]
    ]
end
