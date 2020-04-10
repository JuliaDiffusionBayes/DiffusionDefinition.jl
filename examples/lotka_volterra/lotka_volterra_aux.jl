@diffusion_process LotkaVolterraAux{K,R} begin
    :dimensions
    process --> 2
    wiener --> 2

    :parameters
    (α, β, γ, δ, σ1, σ2) --> K

    :auxiliary_info
    t0 --> Float64
    T --> Float64
    vT --> R

    :additional
    constdiff --> true
    linear --> true
end

DiffusionDefinition.B(t, P::LotkaVolterraAux) = @SMatrix [-0.0 -P.β*P.γ/P.δ; P.α*P.δ/P.β -0.0]
DiffusionDefinition.β(t, P::LotkaVolterraAux) = @SVector [P.γ/P.δ*P.α, -P.α/P.β*P.γ]
DiffusionDefinition.σ(t, P::LotkaVolterraAux) = SDiagonal(P.σ1, P.σ2)
