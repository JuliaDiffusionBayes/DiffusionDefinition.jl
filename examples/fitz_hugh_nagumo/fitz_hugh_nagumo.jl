function DiffusionDefinition.diffusion_message(::Val{:FitzHughNagumo})
    println("* * * * * * * * * * * * *")
    println("FitzHugh–Nagumo diffusion")
    println("* * * * * * * * * * * * *")
    println("Description")
    println("------------")
    println("A model developed to mimic the evolution of a neuron's membrane ")
    println("potential. Represented by two-dimensional process (Yₜ, Xₜ)")
    println("dYₜ = (Yₜ - (Yₜ)³ - Xₜ + s)/ϵ dt,")
    println("dXₜ = (γYₜ - Xₜ - β)dt + σdWₜ.")
    println("The model has five parameters:")
    println("\tϵ, s, γ, β, σ")
    println("\nConstructors")
    println("------------")
    println("You may initialize the process either by passing the parameters")
    println("in the order above, as in:")
    println("ϵ, s, γ, β, σ = ...")
    println("P = FitzHughNagumo(ϵ, s, γ, β, σ)")
    println("or by keyword arguments")
    println("p = (ϵ=..., s=..., γ=..., β=..., σ=...)")
    println("P = FitzHughNagumo(;p...)")
    println()
end

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
