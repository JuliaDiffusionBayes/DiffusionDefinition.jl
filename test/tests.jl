@diffusion_process Prokaryote begin
    :dimensions
    #---------
    process --> 3
    wiener --> 5

    :parameters
    #---------
    _ --> (3, Float64)
    param --> Float64
    stem --> (2, Int64)
    (theta, alpha) --> (2, Int32)
    (beta, gamma) --> Float32
    (yota, zeta) --> (Float64, Int64)
    v --> SArray{Tuple{4},Float64,1,4}

    :conjugate
    #---------
    phi(t, u) --> (
        (p1*u[2]*(u[2]-1),u[1],u[1]),
        (p3*u[5]*t,u[2],u[4]),
    )
    nonhypo(x) --> x
    num_non_hypo --> 4

    :additional
    #---------
    domain --> LowerBoundedDomain((0.0, 1.0, 0.0, 0.0), (1,2,3,4))
    constdiff --> false
    linear --> false
    eltype --> Float64
end a b c d

struct ABC
    x::Int64
end

temp = ABC(1)
Meta.parse("temp.x").args
eval(Expr(:., :temp, QuoteNode(:x)))


print(repr(quote ($(Expr(:., :P, :p1)) * x[5] * t, x[2], x[4]) end))

m = quote $(
    Expr(:-->,
        :(
            phi(::Val{0}, t, x, P::Prokaryote)
        ),
        :(
            ($(Expr(:., :P, :p1)) * x[2] * (x[2] - 1), x[1], x[1])
        )
    )
) end

phi_def = Expr(:call,
    :phi,
    Expr(:(::),
        Val{0}
    ),
    :t,
    :x,
    Expr(:(::), :P, :Prokaryote)
)

phi_body = Expr(:tuple,
    Expr(:call,
        :*,
        Expr(:., :P, :p1),
        :x,
        :t
    ),
    :x,
    :t,
)

eval(
    Expr(:(=), phi_def, phi_body)
)
Expr(:tuple).args

m1 = Meta.parse("""phi(::Val{0}, t, x, P::Prokaryote) = (P.p1 * x[2] * (x[2] - 1), x[1], x[1])""")
Expr(m1.head, m1.args[1], m1.args[2])

m1.head
m1.args[1]
m1.args[2].head
#m1.args[2].args[1]
m1.args[2].args[2].args[1].head

eval(Expr(m1.head,
    m1.args[1],
    Expr(:block,
        #m1.args[2].args[1],
        m1.args[2].args[2],
    )
))

repr(m1)

blocktemp = :(begin end)
push!(blocktemp.args, m1.args[2].args)
eval(Expr(:-->, m1.args[1], blocktemp))
eval(m)

m = Meta.parse("""
phi(P::Prokaryote, x, ::Val{0}) --> (
    P.p1*x[2]*(x[2]-1),
    p3*x[5]*t,
)
""")
m.args[1].args[4].args[1].args

m.args[1].args[2].args
phi_sym = m.args[1].args[1]
t_sym = m.args[1].args[2]
x_sym = m.args[1].args[3]

m.args[2].head
m.args[2].args[1].args[3].args[1] = :u
m
