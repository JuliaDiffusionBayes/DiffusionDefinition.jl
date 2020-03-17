@diffusion_process FitzHughDiffusion begin
    :dimensions
    process --> 2
    wiener --> 1

    :parameters
    (ϵ, s, γ, β, σ) --> Float64

    :additional
    constdiff --> true
end

function b(t, x, P::FitzHughNagumo)
    ℝ{2}(
        (x[1] - x[2] - x[1]^3 + P.s)/P.ϵ,
        P.γ*x[1] - x[2] + P.β
    )
end

σ(t, x, P::FitzHughNagumo) = ℝ{2}(0.0, P.σ)
