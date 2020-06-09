function DiffusionDefinition.diffusion_message(::Val{:LotkaVolterra})
    println("* * * * * * * * * * *")
    println("Lotka–Volterra system")
    println("* * * * * * * * * * *")
    println("Description")
    println("------------")
    println("A simple, scalar-valued predator-prey model.")
    println("It is represented by a two-dimensional process (Xₜ, Yₜ).")
    println("dXₜ = (αXₜ - βXₜYₜ)dt + σ₁dW¹ₜ,")
    println("dYₜ = (δXₜYₜ - γYₜ)dt + σ₂dW²ₜ.")
    println("The model has six parameters:")
    println("\tα, β, γ, δ, σ1, σ2")
    println("\nConstructors")
    println("------------")
    println("You may initialize the process either by passing the parameters")
    println("in the order above, as in:")
    println("α, β, γ, δ, σ1, σ2 = ...")
    println("P = LotkaVolterra(α, β, γ, δ, σ1, σ2)")
    println("or by keyword arguments")
    println("p = (α=..., β=..., γ=..., δ=..., σ1=..., σ2=...)")
    println("P = LotkaVolterra(;p...)")
    println()
end

@diffusion_process LotkaVolterra{K} begin
    :dimensions
    process --> 2
    wiener --> 2

    :parameters
    (α, β, γ, δ, σ1, σ2) --> K

    #=
    :conjugate
    phi(t, x) --> (
        (0.0, 0.0),
        (x[1], 0.0),
        (-x[1]*x[2], 0.0),
        (0.0, x[1]*x[2]),
        (0.0, -x[2]),
        (0.0, 0.0),
        (0.0, 0.0)
    )
    nonhypo(x) --> x
    num_non_hypo --> 2
    =#

    :additional
    constdiff --> true
    statespace --> LowerBoundedStateSpace((1,2),(0.0, 0.0))
end

function DiffusionDefinition.b(t, x, P::LotkaVolterra)
    @SVector [
        P.α*x[1] - P.β*x[1]*x[2],
        P.δ*x[1]*x[2] - P.γ*x[2],
    ]
end

DiffusionDefinition.σ(t, x, P::LotkaVolterra) = SDiagonal(P.σ1, P.σ2)
