@diffusion_process FitzHughNagumoAux{K,R} begin
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

B(t, P::FitzHughNagumoAux) = @SMatrix [
    1/P.ϵ-3*P.vT[1]^2/P.ϵ   -1/P.ϵ;
    P.γ                 -1.0
]

β(t, P::FitzHughNagumoAux) = @SVector [
    P.s/P.ϵ+2*P.vT[1]^3/P.ϵ,
    P.β
]

σ(t, P::FitzHughNagumoAux) = @SVector [0.0, P.σ]
