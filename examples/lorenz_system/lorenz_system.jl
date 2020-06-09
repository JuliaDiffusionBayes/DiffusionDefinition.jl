function DiffusionDefinition.diffusion_message(::Val{:Lorenz})
    println("* * * * * * *")
    println("Lorenz system")
    println("* * * * * * *")
    println("Description")
    println("------------")
    println("A three-dimensional elliptic diffusion")
    println("dXₜ = θ₁(Yₜ-Xₜ)dt + σdW¹ₜ,")
    println("dYₜ = [Xₜ(θ₂ - Zₜ)-Yₜ]dt + σdW²ₜ,")
    println("dZₜ = (XₜYₜ - θ₃Zₜ)dt + σdW³ₜ.")
    println("The model has four parameters:")
    println("\tθ1, θ2, θ3, σ")
    println("\nConstructors")
    println("------------")
    println("You may initialize the process either by passing the parameters")
    println("in the order above, as in:")
    println("θ₁, θ₂, θ₃, σ = ...")
    println("P = Lorenz(θ₁, θ₂, θ₃, σ)")
    println("or by keyword arguments")
    println("p = (θ1=..., θ2=..., θ3=..., σ=...)")
    println("P = Lorenz(;p...)")
    println()
end

@diffusion_process Lorenz{T} begin
    :dimensions
    process --> 3
    wiener --> 3

    :parameters
    _ --> (3, T)
    σ --> T

    :additional
    constdiff --> true
    sparsediff --> true
end

function DiffusionDefinition.b(t, x, P::Lorenz)
    @SVector [
        P.p1*(x[2]-x[1]),
        P.p2*x[1] - x[2] - x[1]*x[3],
        x[1]*x[2] - P.p3*x[3]
    ]
end

DiffusionDefinition.σ(t, x, P::Lorenz) = SDiagonal(P.σ, P.σ, P.σ)
