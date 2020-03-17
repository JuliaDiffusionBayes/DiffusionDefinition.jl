@diffusion_process Sine begin
    :parameters
    (a, b, c, σ) --> Float64

    :additional
    constdiff --> true
end

b(t, x, P::Sine) = ℝ{1}(P.a + P.b*sin.(P.c * x))

σ(t, x, P::Sine) = ℝ{1}(P.σ)
