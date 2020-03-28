@diffusion_process LotkaVolterraAux begin
    :dimensions
    process --> 2
    wiener --> 2

    :parameters
    (α, β, γ, δ, σ1, σ2, t, T) --> Float64
    (u, v) --> ℝ{2}

    :conjugate
    phi(t, x) --> (
        (0.0, 0.0),
        (x[1], 0.0),
        (-x[1]*x[2], 0.0),
        (0.0, x[1]*x[2]),
        (0.0, -x[2]),
        (0.0, 0.0),
        (0.0, 0.0)
    )
    nonhypo(x) --> x
    num_non_hypo --> 2

    :additional
    constdiff --> true
    linear --> true
end

DiffusionDefinition.B(t, P::LotkaVolterraAux) = @SMatrix [-0.0 -P.β*P.γ/P.δ; P.α*P.δ/P.β -0.0]
DiffusionDefinition.β(t, P::LotkaVolterraAux) = ℝ{2}(P.γ/P.δ*P.α, -P.α/P.β*P.γ)
DiffusionDefinition.σ(t, P::LotkaVolterraAux) = SDiagonal(P.σ1, P.σ2)
