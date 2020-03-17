@diffusion_process FitzHughDiffusionAlt begin
    :dimensions
    process --> 2
    wiener --> 1

    :parameters
    (ϵ, s, γ, β, σ) --> Float64

    :additional
    constdiff --> true
end

function b(t, x, P::FitzHughNagumoAlt)
    @SVector [
        x[2],
        -(
            (P.γ-1.0)*x[1]
            + x[1]^3
            + P.ϵ*x[2]
            - P.s
            + P.β
            + (3.0*x[1]^2 - 1.0)*x[2]
        )/P.ϵ
    ]
end

σ(t, x, P::FitzHughNagumoAlt) = @SVector [0.0, P.σ/P.ϵ]
