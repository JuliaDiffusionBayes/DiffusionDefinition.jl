@diffusion_process LotkaVolterraAux begin
    :dimensions
    process --> 2
    wiener --> 2

    :parameters
    (α, β, γ, δ, σ1, σ2, t, T) --> Float64
    (u, v) --> SArray{Tuple{2},Float64,1,2}

    :additional
    constdiff --> true
    statespace --> LowerBoundedStateSpace((1,2),(0.0, 0.0))
end

B(t, P::LotkaVolterraAux) = @SMatrix [-0.0 -P.β*P.γ/P.δ; P.α*P.δ/P.β -0.0]
β(t, P::LotkaVolterraAux) = @SVector [P.γ/P.δ*P.α, -P.α/P.β*P.γ]
σ(t, P::LotkaVolterraAux) = SDiagonal(P.σ1, P.σ2)
