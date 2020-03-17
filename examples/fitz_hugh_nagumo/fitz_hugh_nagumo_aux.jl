@diffusion_process FitzHughDiffusionAux begin
    :dimensions
    process --> 2
    wiener --> 1

    :parameters
    (ϵ, s, γ, β, σ, t, T) --> Float64
    (u, v) --> ℝ{2}

    :additional
    constdiff --> true
    linear --> true
end

function B(t, P::FitzHughDiffusionAux)
    @SMatrix [
        1/P.ϵ-3*P.v^2/P.ϵ   -1/P.ϵ;
        P.γ                 -1.0
    ]
end

function β(t, P::FitzHughDiffusionAux)
    ℝ{2}(
        P.s/P.ϵ+2*P.v^3/P.ϵ,
        P.β
    )
end

σ(t, P::FitzHughDiffusionAux) = ℝ{2}(0.0, P.σ)
