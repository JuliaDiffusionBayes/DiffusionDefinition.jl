function DiffusionDefinition.diffusion_message(::Val{:Sine})
    println("* * * * * * * ")
    println("Sine diffusion")
    println("* * * * * * * ")
    println("Description")
    println("------------")
    println("A one dimensional diffusion Xₜ:")
    println("dXₜ=(a + b×sin(c Xₜ))dt + σdWₜ.")
    println("Xₜ ∈ ℝ.")
    println("The model has four parameters:")
    println("\ta, b, c, σ")
    println()
    println("Constructors")
    println("------------")
    println("You may initialize the process either by passing the parameters")
    println("in the order above, as in:")
    println("a, b, c, σ = ...")
    println("P = Sine(a, b, c, σ)")
    println("or by keyword arguments")
    println("p = (a = ..., b = ..., c = ..., σ = ...)")
    println("P = Sine(;p...)")
    println()
end

@diffusion_process Sine{T} begin
    :parameters
    (a, b, c, σ) --> T

    :additional
    constdiff --> true
end

DiffusionDefinition.b(t, x, P::Sine) = P.a + P.b*sin(P.c * x)
DiffusionDefinition.σ(t, x, P::Sine) = P.σ

DiffusionDefinition.default_type(::Sine) = Float64
DiffusionDefinition.default_wiener_type(::Sine) = Float64
