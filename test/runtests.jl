using DiffusionDefinition
using Test
using StaticArrays, Random
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

        :aux_info
        T --> Float64
        vT --> SArray{Tuple{4},Float64,1,4}

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

    param_tuple = (1.0, 2.0, 3.0, 4.0, 5, 6, Int32(7), Int32(8),
                      Float32(9.0), Float32(10.0), 11.0, 12)
    end_point_info_tuple = (1.0, (@SVector [1.0, 2.0, 3.0, 4.0]))
    P = TestDiffusion(param_tuple..., end_point_info_tuple...)

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

    param_names = (:p1, :p2, :p3, :param, :stem1, :stem2, :theta, :alpha, :beta, :gamma, :yota, :zeta, :T, :vT, :xT)
    @test DD.parameter_names(P) == param_names
    #TODO introduce this back ?
    #@test DD.parameter_names(typeof(P)) == (:p1, :p2, :p3, :param, :stem1, :stem2, :theta, :alpha, :beta, :gamma, :yota, :zeta)
    params_is = DD.parameters(P)
    params_should = Dict([n=>p for (n,p) in zip(param_names, param_tuple)]...)
    for p in keys(params_should)
        @test params_is[p] == params_should[p]
    end
    @test DD.end_point_info_names(P) == (:T, :vT, :xT)
    epi_is = DD.end_point_info(P)
    epi_should = Dict(
        :T=>end_point_info_tuple[1],
        :vT=>end_point_info_tuple[2],
        :xT=>(@SVector [0.0, 0.0, 0.0])
    )
    for p in keys(epi_should)
        @test epi_is[p] == epi_should[p]
    end

    @diffusion_process TestDiffusion2 begin
        :dimensions
        process --> 3
        wiener --> 5

        :parameters
        _ --> (3, Float64)

        :aux_info
        T --> Float64
        vT --> SArray{Tuple{4},Float64,1,4}

    end

    param_tuple = (1.0, 2.0, 3.0)
    end_point_info_tuple = (1.0, (@SVector [1.0, 2.0, 3.0, 4.0]))

    P = TestDiffusion2(param_tuple..., end_point_info_tuple...)

    @test DD.parameter_names(P) == (:p1, :p2, :p3, :T, :vT, :xT)

    params_is = DD.parameters(P)
    params_should = Dict(:p1 => 1.0, :p2 => 2.0, :p3 => 3.0, :T => 1.0, :vT => (@SVector [1.0, 2.0, 3.0, 4.0]))
    for p in keys(params_should)
        @test params_is[p] == params_should[p]
    end
    @test DD.end_point_info_names(P) == (:T, :vT, :xT)
    epi_is = DD.end_point_info(P)
    epi_should = Dict(
        :T=>end_point_info_tuple[1],
        :vT=>end_point_info_tuple[2],
        :xT=>(@SVector [0.0, 0.0, 0.0])
    )
    for p in keys(epi_should)
        @test epi_is[p] == epi_should[p]
    end
    P1 = TestDiffusion2(1.0, 1.0, 1.0, end_point_info_tuple...)
    P2 = TestDiffusion2(1.0, 1.0, 1.0, 1.0, (@SVector [1.0, 1.0, 1.0, 1.0]), (@SVector [1.0, 1.0, 1.0]))

    #@test DD.update_params(P, [1.0, 1.0, 1.0]...) == P1

    #@test DD.clone(P, [1.0, 1.0, 1.0]..., [1.0, (@SVector [1.0, 1.0, 1.0, 1.0]), (@SVector [1.0, 1.0, 1.0])]...) == P2
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

@testset "utility functions" begin
    @test DD.remove_curly(Vector{Float64}) == Array
    @test DD.get_curly(Vector{Float64}) == (:Float64, 1)
end

@testset "Buffers" begin
    d, m = 4,5
    el = Float64
    sb = DD.StandardEulerBuffer(zeros(el,d), zeros(el,(d,m)), zeros(el,m))
    #@test size(sb) == (d + d + d*m + m,)
    @test length(sb.b) == d
    @test size(sb.σ) == (d, m)
    @test length(sb.dW) == m
    #@test eltype(sb) == el
    #sb2 = similar(sb)
    #@test sb2 !== sb
    #sb3 = similar(sb, SVector{3,Float64})
    #@test size(sb3) == (d + d*m + m,)
    #@test sb3.data == zeros(SVector{3,Float64}, d+d*m+m)
    #=
    @testset "indexing StandardEulerBuffer" begin
        sb.b[1] = 1.0
        #@test sb[:] == sb.data
        sb[2:5] .= [2.0, 3.0, 4.0, 5.0]
        sb[[6,8]] .= [6.0, 8.0]
        sb[15] = 15.0
        sb[25] = 25.0

        @test sb.b == [1.0, 2.0, 3.0, 4.0]
        @test sb.y == [5.0, 6.0, 0.0, 8.0]
        @test sb.σ == [
            0.0 0.0 0.0 0.0 25.0;
            0.0 0.0 0.0 0.0 0.0;
            0.0 15.0 0.0 0.0 0.0;
            0.0 0.0 0.0 0.0 0.0
        ]
        @test sb.dW == [0.0, 0.0, 0.0, 0.0, 0.0]
    end
    =#
    lb = DD.LinearDiffBuffer(zeros(el,d), zeros(el,(d,m)), zeros(el,m), zeros(el,(d,d)))
    #@test size(lb) == (d + d + d*m + d*d + m,)
    @test length(lb.b) == d
    @test size(lb.B) == (d,d)
    @test size(lb.σ) == (d, m)
    @test length(lb.dW) == m
end

@testset "setting parameters and cloning" begin
    @load_diffusion LotkaVolterra
    @load_diffusion LotkaVolterraAux

    θ = [2.0/3.0, 4.0/3.0, 1.0, 1.0, 0.2, 0.3]

    P_target = LotkaVolterra(θ...)
    P_aux = LotkaVolterraAux(θ..., 0.0, 1.0, (@SVector [1.0, 2.0]))

    p_names = (:α, :β, :γ, :δ, :σ1, :σ2)
    pars = Dict([pn=>p for (pn,p) in zip(p_names,θ)])
    ep_names = (:t0, :T, :vT, :xT)

    @test DD.parameter_names(LotkaVolterra) == p_names
    @test DD.parameter_names(P_target) == p_names
    @test DD.parameters(P_target) == pars

    @test DD.const_parameter_names(LotkaVolterra) == tuple()
    @test DD.const_parameter_names(P_target) == tuple()
    @test DD.const_parameters(P_target) == Dict{Any,Any}()

    @test DD.var_parameter_names(LotkaVolterra) == p_names
    @test DD.var_parameter_names(P_target) == p_names
    @test DD.var_parameters(P_target) == pars

    @test DD.end_point_info_names(LotkaVolterra) == tuple()
    @test DD.end_point_info_names(P_target) == tuple()
    @test DD.end_point_info(P_target) == Dict{Any,Any}()

    @test DD.parameter_names(LotkaVolterraAux) == Tuple(vcat(p_names..., ep_names...))
    @test DD.parameter_names(P_aux) == Tuple(vcat(p_names..., ep_names...))

    #TODO finish below...
    DD.parameters(P_aux)

    DD.const_parameter_names(LotkaVolterraAux)
    DD.const_parameter_names(P_aux)
    DD.const_parameters(P_aux)

    DD.var_parameter_names(LotkaVolterraAux)
    DD.var_parameter_names(P_aux)
    DD.var_parameters(P_aux)

    DD.end_point_info_names(LotkaVolterraAux)
    DD.end_point_info_names(P_aux)
    DD.end_point_info(P_aux)

    DD.clone(
        P_target,
        [10.0, 20.0, 30.0, 40.0],
        Dict(2=>2, 4=>4),
        [
            (global_idx=2, pname=:δ),
            (global_idx=4, pname=:α)
        ],
        Val(:associate_by_name)
    )
end


@testset "Sampling trajectories" begin

    Random.seed!(3)
    N, dt, T = 30, 0.0001, 10.0
    tt = collect(0.0:dt:T)
    quad_var(x) = sum(abs2.(diff(x)))

    tr = trajectory(tt, SVector{N,Float64})
    wr = Wiener()
    rand!(wr, tr)

    tr_mut = trajectory(tt, Vector{Float64}, N)
    wr_mut = Wiener(N)
    rand!(wr_mut, tr_mut)

    @test maximum([abs(quad_var(map(x->x[i], tr.x))/T-1.0) for i in 1:N])<0.05
    @test maximum([abs(quad_var(map(x->x[i], tr_mut.x))/T-1.0) for i in 1:N])<0.05

    # an example of a pre-defined diffusion
    @load_diffusion LotkaVolterraAux

    α, β, γ, δ, σ1, σ2 = 2.0/3.0, 4.0/3.0, 1.0, 1.0, 0.2, 0.3
    lv_aux = LotkaVolterraAux(
        α, β, γ, δ, σ1, σ2,
        0.0, 1.0, zero(DD.ℝ{2}), zero(DD.ℝ{2})
    )
    N, dt, T = 2, 0.0001, 1000.0
    tt = collect(0.0:dt:T)

    WW = trajectory(tt, SVector{N,Float64})
    wr = Wiener()
    rand!(wr, WW)

    v0 = @SVector [1.0, 2.0]
    XX = trajectory(tt, SVector{N,Float64})
    DD.solve!(XX, WW, lv_aux, v0)

    mean(x) = sum(x)/length(x)
    x1 = map(x->x[1], XX.x)
    x2 = map(x->x[2], XX.x)

    @test abs(mean(x1)-γ/δ) < 0.05
    @test abs(mean(x2)-α/β) < 0.05
    @test abs(quad_var(x1)/T-σ1^2) < 0.01
    @test abs(quad_var(x2)/T-σ2^2) < 0.01

    # and in-place simulation
    function DD.B!(buffer, t, P::LotkaVolterraAux)
        #B = DD.get_B(buffer)
        buffer.B[1,1] = 0.0
        buffer.B[1,2] = -P.β*P.γ/P.δ
        buffer.B[2,1] = P.α*P.δ/P.β
        buffer.B[2,2] = 0.0
    end

    function DD.β!(buffer, t, P::LotkaVolterraAux)
        #b = DD.get_b(buffer)
        buffer.b[1] = P.γ/P.δ*P.α
        buffer.b[2] = -P.α/P.β*P.γ
    end

    function DD.σ!(buffer, t, P::LotkaVolterraAux)
        #σ = DD.get_σ(buffer)
        buffer.σ[1,1] = P.σ1
        buffer.σ[1,2] = 0.0
        buffer.σ[2,1] = 0.0
        buffer.σ[2,2] = P.σ2
    end

    WW = trajectory(tt, Vector{Float64}, N)
    wr = Wiener(N)
    rand!(wr, WW)

    XX = trajectory(tt, Vector{Float64}, N)
    v0 = [1.0, 2.0]
    buffer = DD.LinearDiffBuffer{Vector{Float64}}(lv_aux)

    DD.solve!(XX, WW, lv_aux, v0, buffer)

    x1 = map(x->x[1], XX.x)
    x2 = map(x->x[2], XX.x)
    @test abs(mean(x1)-γ/δ) < 0.05
    @test abs(mean(x2)-α/β) < 0.05
    @test abs(quad_var(x1)/T-σ1^2) < 0.05
    @test abs(quad_var(x2)/T-σ2^2) < 0.05
end
