using DiffusionDefinition
using Test
using StaticArrays
DD = DiffusionDefinition

@testset "DiffusionDefinition.jl" begin
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
        statespace --> LowerBoundedStateSpace((1,2,3), (0.0, 1.0, 0.0))
        constdiff --> false
        linear --> false
        eltype --> Float64
    end a b c d


    P = TestDiffusion(1.0, 2.0, 3.0, 4.0, 5, 6, Int32(7), Int32(8),
                      Float32(9.0), Float32(10.0), 11.0, 12,
                      @SVector [1.0, 2.0, 3.0, 4.0])

    @test DD.nonhypo(P, :x) == :x
    @testset "satisfying bounds" begin
        @test !DD.bound_satisfied(P, [0.1, 2.0, 0.0])
        @test DD.bound_satisfied(P, [0.1, 2.0, 0.1])
        @test !DD.bound_satisfied(P, [0.1, 1.0, 0.1])
    end
    @test !DD.constdiff(P)
    @test typeof(P) <: DD.DiffusionProcess
    @test !(typeof(P) <: DD.LinearDiffusion)
    @test eltype(P) == Float64
end


@testset "State space restrictions" begin
    lb = LowerBoundedStateSpace((1,3), (-2.0, 3.0))
    ub = UpperBoundedStateSpace((1,3), (4.0, 10.0))
    b = BoundedStateSpace(lb, ub)
    nb = UnboundedStateSpace()

    a = [-1,0,1.0,4.0]
    v = [-1,0,1.0,20.0]
    @test DD.bound_satisfied(nb, a)
    @test !DD.bound_satisfied(lb, a)
    @test DD.bound_satisfied(ub, a)
    @test !DD.bound_satisfied(b, a)
    @test DD.bound_satisfied(nb, v)
    @test !DD.bound_satisfied(lb, v)
    @test DD.bound_satisfied(ub, v)
    @test !DD.bound_satisfied(b, v)
end
