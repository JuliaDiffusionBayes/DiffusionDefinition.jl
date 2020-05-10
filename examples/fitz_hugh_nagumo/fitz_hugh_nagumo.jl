@diffusion_process FitzHughNagumo{T} begin
    :dimensions
    process --> 2
    wiener --> 1

    :parameters
    (ϵ, s, γ, β, σ) --> T

    :additional
    constdiff --> true
end

DiffusionDefinition.b(t, x, P::FitzHughNagumo) = @SVector [
    (x[1] - x[2] - x[1]^3 + P.s)/P.ϵ,
    P.γ*x[1] - x[2] + P.β
]

DiffusionDefinition.σ(t, x, P::FitzHughNagumo) = @SMatrix [0.0; P.σ]
