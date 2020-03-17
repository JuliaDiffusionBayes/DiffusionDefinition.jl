@diffusion_process LotkaVolterra begin
    :dimensions
    process --> 2
    wiener --> 2

    :parameters
    (α, β, γ, δ, σ1, σ2) --> Float64

    :additional
    constdiff --> true
    statespace --> LowerBoundedStateSpace((1,2),(0.0, 0.0))
end

function b(t, x, P::LotkaVolterra)
    ℝ{2}(
        P.α*x[1] - P.β*x[1]*x[2],
        P.δ*x[1]*x[2] - P.γ*x[2],
    )
end

σ(t, x, P::LotkaVolterra) = SDiagonal(P.σ1, P.σ2)
