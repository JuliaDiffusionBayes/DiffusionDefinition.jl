@diffusion_process FitzHughDiffusionAltAuxLin begin
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

function B(t, P::FitzHughDiffusionAltAuxLin)
    @SMatrix [
        0.0  1.0;
        (1.0-P.γ-3.0*P.v[1]^2-6*P.v[1]*P.v[2])/P.ϵ (1.0-P.ϵ-3.0*P.v[1]^2)/P.ϵ
    ]
end

function β(t, P::FitzHughDiffusionAltAuxLin)
    ℝ{2}(
        0.0,
        (2*P.v[1]^3+P.s-P.β+6*P.v[1]^2*P.v[2])/P.ϵ
    )
end

σ(t, P::FitzHughDiffusionAltAuxLin) = ℝ{2}(0.0, P.σ/P.ϵ)
