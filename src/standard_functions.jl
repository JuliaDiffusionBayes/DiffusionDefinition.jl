"""
    dimension(::DiffusionProcess{T,DP,DW})

Return dimension of the stochastic process and driving Brownian motion.
"""
dimension(::DiffusionProcess{T,DP,DW}) where {T,DP,DW} = (
    process = DP,
    wiener = DW
)

dim_process(::DiffusionProcess{T,DP}) where {T,DP} = DP
dim_wiener(::DiffusionProcess{T,DP,DW}) where {T,DP,DW} = DW

"""
    Base.eltype(::DiffusionProcess{T}) where T = T

Return the datatype that each coordinate of the stochastic process is stored in.
"""
Base.eltype(::DiffusionProcess{T}) where T = T

"""
    state_space(::DiffusionProcess{T,DP,DW,SS})

Return the state space restrictions.
"""
state_space(::DiffusionProcess{T,DP,DW,SS}) where {T,DP,DW,SS} = SS


"""
    default_type(::DiffusionProcess{T,DP})

Allows for inference of data type that encodes the state space of a given
diffusion.
"""
default_type(::DiffusionProcess{T,DP}) where {T,DP} = SVector{DP,T}

"""
    default_wiener_type(::DiffusionProcess{T,DP,DW})

Allows for inference of data type that encodes the state space of the Brownian
motion driving a given diffusion process.
"""
default_wiener_type(::DiffusionProcess{T,DP,DW}) where {T,DP,DW} = SVector{DW,T}

"""
    parameters(P::DiffusionProcess)

Return a tuple of pairs of `parameter_name` => `parameter_value`.
"""
function parameters(P::DiffusionProcess) end

"""
    const_parameters(P::DiffusionProcess)

Return a tuple of pairs of `parameter_name` => `parameter_value`. Return only
those parameteres that are considered to be `constant`.
"""
function const_parameters(P::DiffusionProcess) end

"""
    var_parameters(P::DiffusionProcess)

Return a tuple of pairs of `parameter_name` => `parameter_value`. Return only
those parameteres that are considered to be `variable`.
"""
function var_parameters(P::DiffusionProcess) end

"""
    parameter_names(P::DiffusionProcess)

Return a tuple with the names of all paremeters.
"""
parameter_names(P::DiffusionProcess) = first.(parameters(P))

"""
    const_parameter_names(P::DiffusionProcess)

Return a tuple with the names of all paremeters that are considered to be
`constant`.
"""
const_parameter_names(P::DiffusionProcess) = first.(const_parameters(P))

"""
    var_parameter_names(P::DiffusionProcess)

Return a tuple with the names of all paremeters that are considered to be
`variable`.
"""
var_parameter_names(P::DiffusionProcess) =first.(var_parameters(P))

a(t, x, P::DiffusionProcess) = σ(t, x, P) * σ(t, x, P)'


#------------------------------------------------------------------------------#
#                     Default fallbacks for linear diffusions
#------------------------------------------------------------------------------#
function b!(buffer, t, x, P::LinearDiffusion)
    B!(buffer, t, P)
    β!(buffer, t, P)
    mul!(buffer.b, buffer.B, x, true, true)
end

σ!(buffer, t, x, P::LinearDiffusion) = σ!(buffer, t, P)
a!(buffer, t, x, P::LinearDiffusion) = a!(buffer, t, P)

b(t, x, P::LinearDiffusion) = B(t, P)*x + β(t, P)
σ(t, x, P::LinearDiffusion) = σ(t, P)
a(t, x, P::LinearDiffusion) = a(t, P)
a(t, P::LinearDiffusion) = σ(t, P) * σ(t, P)'

#------------------------------------------------------------------------------#
#
#                   Default fallbacks for indexed time
#
#------------------------------------------------------------------------------#
_b((t,i)::IndexedTime, x, P::DiffusionProcess) = b(t, x, P)
_σ((t,i)::IndexedTime, x, P::DiffusionProcess) = σ(t, x, P)
_b!(buffer, (t,i)::IndexedTime, x, P::DiffusionProcess) = b!(buffer, t, x, P)
_σ!(buffer, (t,i)::IndexedTime, x, P::DiffusionProcess) = σ!(buffer, t, x, P)

"""
    diagonaldiff(P::DiffusionProcess)

Indicator for whether the volatility coefficient is represented by a diagonal
matrix
"""
diagonaldiff(P::DiffusionProcess) = false

"""
    sparsediff(P::DiffusionProcess)

Indicator for whether the volatility coefficient is represented by a sparse
matrix
"""
sparsediff(P::DiffusionProcess) = false

"""
    diagonalBmat(P::DiffusionProcess)

Indicator for whether the B matrix (if exists) is represented by a diagonal
matrix
"""
diagonalBmat(P::DiffusionProcess) = false

"""
    sparseBmat(P::DiffusionProcess)

Indicator for whether the B matrix (if exists) is represented by a sparse
matrix
"""
sparseBmat(P::DiffusionProcess) = false

"""
    clone(P::T, args...)

Instantiate a new diffusion process of the same type as `P` but with new
arguments `args`.
"""
clone(P::T, args...) where T<:DiffusionProcess = T(args...)

function update_params(P::T, new_params...) where T<:DiffusionProcess
    ei = end_point_info(P)
    ei === nothing ? T(new_params...) : T(new_params..., ei...)
end


Base.zero(P::DiffusionProcess) = zero(P, Val(:process))

"""
    Base.zero(P::DiffusionProcess)

Instantiate a zero element that can represent a state of a diffusion.
"""
function Base.zero(P::DiffusionProcess, ::Val{:wiener})
    deft = default_wiener_type(P)
    zero(deft, dim_wiener(P), ismutable(deft))
end

"""
    Base.zero(P::DiffusionProcess)

Instantiate a zero element that can represent a state of a Brownian motion
driving a diffusion process.
"""
function Base.zero(P::DiffusionProcess, ::Val{:process})
    deft = default_type(P)
    zero(deft, dim_process(P), ismutable(deft))
end

"""
    Base.zero(K::Type, D, ::Val)

If `K` is a mutable type, then create `zeros` of dimension `D` and entries with
types `eltype(K)`. Otherwise, calls regular zero(`K`).
"""
Base.zero

Base.zero(K::Type, D, ::Val{true}) = zeros(eltype(K), D)
Base.zero(K::Type, D, ::Val{false}) = zero(K)



"""
    end_point_info

Return information about the end-point (works only if some information of this
kind has been passed at the time of defining a struct) TODO improve
"""
function end_point_info end

"""
    end_point_info_names

Return names of information pieces about the end-point (works only if some
information of this kind has been passed at the time of defining a struct)
TODO improve
"""
function end_point_info_names end


#custom_zero(D::Integer, ::Type{K}) where K <: Array = zeros(eltype(K), D)
#custom_zero(D::Integer, ::Type{K}) where K = zero(K)

#=
function sequential_add!(save_to, x, y, a=true, b=true, c=false)
    for i in eachindex(save_to, x, y)
        @inbounds save_to[i] = a*x[i]+b*y[i]+c
    end
    return nothing
end

function sequential_mul!(save_to, x, y, a=true, b=true)
    for i in eachindex(save_to, x, y)
        @inbounds save_to[i] = a*x[i]*y[i]+b
    end
    return nothing
end

function sequential_mul_add!(save_to, X, Y, Z=Y, a=true, b=false, c=false)
    for i in eachindex(save_to, X, Y, Z)
        @inbounds save_to[i] = a * X[i] * Y[i] + b * Z[i] + c
    end
    return nothing
end
=#