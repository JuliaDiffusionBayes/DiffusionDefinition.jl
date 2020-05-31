@diffusion_process FitzHughNagumoConjug{T} begin
    :dimensions
    process --> 2
    wiener --> 1

    :parameters
    (ϵ, s, γ, β, σ) --> T

    :additional
    constdiff --> true
end

@conjugate_gaussian FitzHughNagumoConjug begin
    :intercept --> (-x[2],)
    :ϵ --> (x[1]-x[1]^3+(1-3*x[1]^2)*x[2],)
    :s --> (one(x[1]),)
    :γ --> (-x[1],)
    :β --> (-one(x[1]),)
    :hypo_a_inv --> 1.0/P.σ^2
    :nonhypo --> 2:2
end

DiffusionDefinition.b(t, x, P::FitzHughNagumoConjug) = @SVector [
    x[2],
    (
        (P.ϵ - P.γ)*x[1]
        - P.ϵ*(x[1]^3 + (3.0*x[1]^2 - 1.0)*x[2])
        + P.s
        - P.β
        - x[2]
    )
]

σ(t, x, P::FitzHughNagumoConjug) = @SMatrix [0.0; P.σ]
