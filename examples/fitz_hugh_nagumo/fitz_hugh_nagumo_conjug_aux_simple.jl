@diffusion_process FitzHughDiffusionConjugAuxSimple begin
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

B(t, P::FitzHughDiffusionConjugAuxSimple) = @SMatrix [0.0  1.0; 0.0 0.0]

β(t, P::FitzHughDiffusionConjugAuxSimple) = ℝ{2}(0.0, 0.0)

σ(t, P::FitzHughDiffusionAux) = ℝ{2}(0.0, P.σ)
