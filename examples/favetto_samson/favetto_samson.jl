@diffusion_process FavettoSamson{T,K} begin
    :dimensions
    process --> 2
    wiener --> 2

    :parameters
    (α, β, λ, μ, σ1, σ2) --> T
    dose --> K

    :additional
    diagdiff --> true
    constdiff --> true
end

DiffusionDefinition.b(t, x, P::FavettoSamson) = @SVector [
    P.α * P.dose(t) - (P.α +P.β)*x[1] + P.μ*x[2],
    P.λ*x[1] - P.μ*x[2]
]

DiffusionDefinition.σ(t, x, P::FavettoSamson) = SDiagonal(P.σ1, P.σ2)
