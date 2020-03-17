@diffusion_process SineAux begin
    :parameters
    (a, b, c, σ, t, T, u, v) --> Float64

    :additional
    constdiff --> true
end

B(t, P::SinDiffusionAux) = @SMatrix [0 + (t/P.T)/5]
β(t, P::SinDiffusionAux) = ℝ{1}((P.v-P.u)/P.T)

σ(t, x, P::Sine) = ℝ{1}(P.σ)
