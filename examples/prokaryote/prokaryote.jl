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
