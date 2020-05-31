@diffusion_process FitzHughNagumoConjugAuxSimple{K,R} begin
    :dimensions
    process --> 2
    wiener --> 1

    :parameters
    #(ϵ, s, γ, β, σ) --> K
    σ --> K

    :auxiliary_info
    t0 --> Float64
    T --> Float64
    vT --> R

    :additional
    constdiff --> true
    linear --> true
end

DiffusionDefinition.B(t, P::FitzHughNagumoConjugAuxSimple) = @SMatrix [
    0.0  1.0;
    0.0  0.0
]

DiffusionDefinition.β(t, P::FitzHughNagumoConjugAuxSimple) = @SVector [
    0.0, 0.0
]

DiffusionDefinition.σ(t, P::FitzHughNagumoConjugAuxSimple) = @SMatrix [
    0.0;
    P.σ
]
