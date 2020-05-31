@diffusion_process FitzHughNagumoConjugAuxLin{K,R} begin
    :dimensions
    process --> 2
    wiener --> 1

    :parameters
    (ϵ, s, γ, β, σ) --> K

    :auxiliary_info
    t0 --> Float64
    T --> Float64
    vT --> R

    :additional
    constdiff --> true
    linear --> true
end

DiffusionDefinition.B(t, P::FitzHughNagumoConjugAuxLin) = @SMatrix [
    0.0  1.0;
    (P.ϵ-P.γ-3.0*P.ϵ*P.v[1]^2-6*P.ϵ*P.v[1]*P.v[2]) (P.ϵ-1.0-3.0*P.ϵ*P.v[1]^2)
]

DiffusionDefinition.β(t, P::FitzHughNagumoConjugAuxLin) = @SVector [
    0.0,
    2*P.ϵ*P.v[1]^3+P.s-P.β+6*P.ϵ*P.v[1]^2*P.v[2]
]

DiffusionDefinition.σ(t, P::FitzHughNagumoConjugAuxLin) = @SMatrix [0.0; P.σ]
