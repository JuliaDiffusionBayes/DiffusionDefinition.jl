function DiffusionDefinition.diffusion_message(::Val{:SIR})
    println("* * * * * * * * * * * * * * * * * * * * * * * * ")
    println("Susceptible-Infectious-Recovered (SIR) diffusion")
    println("* * * * * * * * * * * * * * * * * * * * * * * * ")
    println("Description")
    println("------------")
    println("A two dimensional diffusion (Iₜ,Rₜ):")
    println("dIₜ=[α(1-Iₜ-Rₜ)Iₜ - βIₜ]dt -σ₁√[(1-Iₜ-Rₜ)Iₜ]dWₜ,")
    println("dIₜ=βIₜdt -σ₂√IₜdWₜ,")
    println("corresponding to proportion of Infectious (Iₜ) and Recovered (Rₜ).")
    println("a proportion of susceptible is computed via: Sₜ=1-Iₜ-Rₜ.")
    println("Iₜ,Rₜ,Sₜ∈[0,1].")
    println("The model has four parameters:")
    println("\tα, β, σ1, σ2")
    println()
    println("Constructors")
    println("------------")
    println("You may initialize the process either by passing the parameters")
    println("in the order above, as in:")
    println("α, β, σ₁, σ₂ = ...")
    println("P = SIR(α, β, σ₁, σ₂)")
    println("or by keyword arguments")
    println("p = (α = ..., β = ..., σ1 = ..., σ2 = ...)")
    println("P = SIR(;p...)")
    println()
end

@diffusion_process SIR{T} begin
    :dimensions
    process --> 2
    wiener --> 2

    :parameters
    (α, β, σ1, σ2) --> T

    #=
    :conjugate
    phi(t, u) --> (
        (0.0, 0.0),
        ((1.0 - u[1] - u[2])*u[1], 0.0),
        (-u[1], u[1]),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
    )
    nonhypo(x) --> x
    num_non_hypo --> 2
    =#

    :additional
    statespace --> BoundedStateSpace(((1,2), (0.0,0.0)), ((1,2), (1.0,1.0)))
end

function DiffusionDefinition.b(t, x, P::SIR)
    @SVector[
        P.α*(1.0 - x[1] - x[2])*x[1] - P.β*x[1],
        P.β*x[1]
    ]
end

function DiffusionDefinition.σ(t, x, P::SIR)
    @SMatrix [
        -P.σ1*sqrt((1.0 - x[1] - x[2])*x[1])  -P.σ2*sqrt(x[1]);
        0.0   P.σ2*sqrt(x[1])
    ]
end

function DiffusionDefinition.bound_satisfied(P::SIR, x)
    (x[1] + x[2]) < 1.0 && DiffusionDefinition._bound_satisfied(P, x)
end
