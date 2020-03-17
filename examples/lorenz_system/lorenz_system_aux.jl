@diffusion_process LorenzAux begin
    :dimensions
    process --> 3
    wiener --> 3

    :parameters
    _ --> (3, Float64)
    (t, T) --> Float64
    (u, v) --> ℝ{3}
    σ --> Float64

    :additional
    constdiff --> true
end

B(t, P::LorenzAux) = @SMatrix [
    -P.p1         P.p1        0.0;
     P.p2-P.v[3]  -1.0    -P.v[1];
     P.v[2]      P.v[1]    -P.p3
]

β(t, P::LorenzAux) = ℝ{3}(
    0.0,
    P.v[1]*P.v[3],
    -P.v[1]*P.v[3]
)

σ(t, P::LorenzAux) = SDiagonal(P.σ, P.σ, P.σ)
