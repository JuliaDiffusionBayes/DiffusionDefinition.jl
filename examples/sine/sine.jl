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
