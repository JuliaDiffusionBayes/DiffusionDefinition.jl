function DiffusionDefinition.diffusion_message(::Val{:Prokaryote})
    println("* * * * * * * * * * ")
    println("Prokaryote diffusion")
    println("* * * * * * * * * * ")
    println("Description")
    println("------------")
    println("Chemical Langevin equation for a simple system describing ")
    println("production of a protein that is repressing its own production.")
    println("It is represented by a 4-dimensional diffusion driven by an ")
    println("8-dimensional Wiener process:")
    println("\ndXₜ=S[θ .× h(Xₜ)]dt + S ⊙ √. ( θ .× h(X_t) ) dWₜ,\n")
    println("where the custom operation ⊙:ℝᵈˣᵏ→ℝᵈˣᵏ is defined via:")
    println("\n(M ⊙ μ)ᵢⱼ = Mᵢⱼμⱼ, i=1,…,d; j=1,…,k,\n")
    println("S is the stoichiometry matrix:")
    println("""S=[
    0   0   1   0   0   0  -1   0
    0   0   0   1  -2   2   0  -1
   -1   1   0   0   1  -1   0   0
   -1   1   0   0   0   0   0   0
]"""
)
    println("and the function h is given by:")
    println("\nh(x) = (x₃x₄, K-x₄, x₄, x₁, x₂(x₂-1)/2, x₃, x₁, x₂)ᵀ.\n")
    println("The underlying process X represents:")
    println("X = (RNA, P, P₂, DNA)")
    println("The model has eight parameters:")
    println("\tc1 c2 c3 c4 c5 c6 c7 c8 K")
    println("The first eight parameters represent the rate of chemical reactions")
    println("The last represent the total number of the DNA strands")
    println("\nConstructors")
    println("------------")
    println("You may initialize the process either by passing the parameters")
    println("in the order above, as in:")
    println("c₁, c₂, c₃, c₄, c₅, c₆, c₇, c₈, K = ...")
    println("P = Prokaryote(c₁, c₂, c₃, c₄, c₅, c₆, c₇, c₈, K)")
    println("or by keyword arguments")
    println("p = (c1=..., c2=..., c3=..., c4=..., c5=..., c6=..., c7=..., c8=..., K=...)")
    println("P = Prokaryote(;p...)")
    println()
end

@diffusion_process Prokaryote{T} begin
    :dimensions
    process --> 4
    wiener --> 8

    :parameters
    c --> (8, T)
    K --> Float64

    :additional
    statespace --> LowerBoundedStateSpace((1,2,3,4), (0.0,1.0,0.0,0.0))
end

function bound_satisfied(P::Prokaryote, x)
    DiffusionDefinition._bound_satisfied(P, x) && (P.K > x[4])
end

@inline _aux₁(x, P::Prokaryote) = P.c5*x[2]*(x[2]-1)
@inline _aux₂(x, P::Prokaryote) = P.c6*x[3]
@inline _aux₃(x, P::Prokaryote) = P.c2*(P.K - x[4])
@inline _aux₄(x, P::Prokaryote) = P.c1*x[3]*x[4]
@inline _aux₁₂(x, P::Prokaryote) = 2.0*_aux₂(x, P) - _aux₁(x, P)
@inline _aux₃₄(x, P::Prokaryote) = _aux₃(x, P) - _aux₄(x, P)


function DiffusionDefinition.b(t, x, P::Prokaryote)
    k1 = _aux₁₂(x, P)
    k2 = _aux₃₄(x, P)
    @SVector [
        P.c3*x[4] - P.c7*x[1],
        P.c4*x[1] + k1 - P.c8*x[2],
        k2 - 0.5 * k1,
        k2
    ]
end

function _σ_prokaryote(x, P)
    k₁ = sqrt(0.5*_aux₁(x, P))
    k₂ = sqrt(_aux₂(x, P))
    k₃ = sqrt(_aux₃(x, P))
    k₄ = sqrt(_aux₄(x, P))
    _O = zero(x[1]) # check if this breaks things for gradients

    @SMatrix [
        _O  _O  sqrt(P.c3*x[4]) _O              _O      _O      -sqrt(P.c7*x[1])    _O;
        _O  _O  _O              sqrt(P.c4*x[1]) -2.0*k₁ 2.0*k₂  _O                  -sqrt(P.c8*x[2]);
        -k₄ k₃  _O              _O              k₁      -k₂     _O                  _O;
        -k₄ k₃  _O              _O              _O      _O      _O                  _O
    ]
end

function _a_prokaryote(x, P)
    k₁ = _aux₁(x, P) + 2.0*_aux₂(x, P)
    k₂ = _aux₃(x, P) + _aux₄(x, P)
    s₁ = P.c3*x[4]+P.c7*x[1]
    s₂ = P.c4*x[1]+P.c8*x[2]
    _O = zero(x[1]) # check if this breaks things for gradients

    @SMatrix [
        s₁  _O      _O          _O;
        _O  s₂+2*k₁ -k₁         _O;
        _O  -k₁     0.5*k₁+k₂   k₂;
        _O  _O      k₂          k₂
    ]
end

DiffusionDefinition.σ(t, x, P::Prokaryote) = _σ_prokaryote(x, P)
DiffusionDefinition.a(t, x, P::Prokaryote) = _a_prokaryote(x, P)
