using Revise
using DiffusionDefinition, StaticArrays
DD = DiffusionDefinition


@diffusion_process TestDiffusion begin
    :dimensions
    process --> 3
    wiener --> 5

    :parameters
    _ --> (3, Float64)
    param --> Float64
    stem --> (2, Int64)
    (theta, alpha) --> (2, Int32)
    (beta, gamma) --> Float32
    (yota, zeta) --> (Float64, Int64)
    v --> SArray{Tuple{4},Float64,1,4}

    :conjugate
    phi(t, u) --> (
        (p1*u[2]*(u[2]-1),u[1],u[1]),
        (p3*u[5]*t,u[2],u[4]),
    )
    nonhypo(x) --> x
    num_non_hypo --> 3

    :additional
    domain --> LowerBoundedStateSpace((1,2,3), (0.0, 1.0, 0.0))
    constdiff --> false
    linear --> false
    eltype --> Float64
end a b c d


P = TestDiffusion(1.0, 2.0, 3.0, 4.0, 5, 6, Int32(7), Int32(8),
                  Float32(9.0), Float32(10.0), 11.0, 12,
                  @SVector [1.0, 2.0, 3.0, 4.0])
v = [0.1, 2.0, 0.0]

DD.constdiff(P)
eltype(P)
