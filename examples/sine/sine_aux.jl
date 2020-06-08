@diffusion_process SineAux{K,R,S} begin
    :parameters
    (a, b, c, σ) --> K

    :auxiliary_info
    t0 --> Float64
    T --> Float64
    v0 --> R
    vT --> S

    :additional
    constdiff --> true
    linear --> true
end

DiffusionDefinition.B(t, P::SineAux) = 0.0 + (t/P.T)/5.0
DiffusionDefinition.β(t, P::SineAux) = (P.vT-P.v0)/P.T
DiffusionDefinition.σ(t, x, P::SineAux) = P.σ

DiffusionDefinition.default_type(::SineAux) = Float64
DiffusionDefinition.default_wiener_type(::SineAux) = Float64
