#==============================================================================#
#
#          Structs defining buffers used for in-place computations
#
#==============================================================================#

"""
    AbstractBuffer{T} <: AbstractArray{T,1}

Types inheriting from `AbstractBuffer` define buffers used for various in-place
computations.
"""
abstract type AbstractBuffer end

"""
    struct StandardEulerBuffer{Tb,Tσ,Tdw} <: AbstractBuffer
        b::Tb
        y::Tb
        σ::Tσ
        dW::Tdw
    end

Standard buffer for the Euler-Maruyama simulations. The intermediary data is
stored in `b`, `σ`, `y` and `dW`.
"""
struct StandardEulerBuffer{Tb,Tσ,Tdw} <: AbstractBuffer
    b::Tb
    y::Tb
    σ::Tσ
    dW::Tdw
end

function StandardEulerBuffer(b::Tb, σ::Tσ, dW::Tdw) where {Tb, Tσ, Tdw}
    y = deepcopy(b)
    StandardEulerBuffer{Tb, Tσ, Tdw}(b, y, σ, dW)
end

function StandardEulerBuffer{K}(P::DiffusionProcess{T,DP,DW}) where {K,T,DP,DW}
    # dW always falls on defaults or not?
    dW = zero(K, DW, ismutable(K)) # zero(P, Val(:wiener))
    b = zero(K, DP, ismutable(K))
    σ = _init_mat_for_buffer(
        Val(diagonaldiff(P)),
        Val(sparsediff(P)),
        K,
        tuple(dimension(P)...)
    )
    StandardEulerBuffer(b, σ, dW)
end

function _init_mat_for_buffer(
        ::Val{false}, # diagonal
        ::Val{false}, # sparse
        ::Type{K}, # eltype
        dims
    ) where K
    zeros(eltype(K), dims)
end

function _init_mat_for_buffer(
        ::Val{true}, # diagonal
        ::Any, # sparse
        ::Type{K}, # eltype
        dims
    ) where K
    Diagonal{eltype(K)}(zeros(eltype(K), dims[1]))
end

function _init_mat_for_buffer(
        ::Val{false}, # diagonal
        ::Val{true}, # sparse
        ::Type{K}, # eltype
        dims
    ) where K
    spzeros(eltype(K), dims...)
end

"""
    struct LinearDiffBuffer{Tb,Tσ,Tdw,TB} <: AbstractBuffer
        b::Tb
        y::Tb
        σ::Tσ
        dW::Tdw
        B::TB
    end

A buffer for Euler-Maruyama simulations of linear diffusions. Almost the same
as `StandardEulerBuffer`, but contains additional space for an intermediate
construction of a matrix `B`.
"""
struct LinearDiffBuffer{Tb,Tσ,Tdw,TB} <: AbstractBuffer
    b::Tb
    y::Tb
    σ::Tσ
    dW::Tdw
    B::TB
end

function LinearDiffBuffer(b::Tb, σ::Tσ, dW::Tdw, B::TB) where {Tb, Tσ, Tdw, TB}
    y = deepcopy(b)
    LinearDiffBuffer{Tb, Tσ, Tdw, TB}(b, y, σ, dW, B)
end

function LinearDiffBuffer{K}(P::LinearDiffusion{T,DP,DW}) where {K,T,DP,DW}
    dW = zero(K, DW, ismutable(K))
    b = zero(K, DP, ismutable(K))
    σ = _init_mat_for_buffer(
        Val(diagonaldiff(P)),
        Val(sparsediff(P)),
        K,
        tuple(dimension(P)...)
    )
    B = _init_mat_for_buffer(
        Val(diagonalBmat(P)),
        Val(sparseBmat(P)),
        K,
        (dim_process(P), dim_process(P))
    )
    LinearDiffBuffer(b, σ, dW, B)
end


"""
    remove_curly(::Type{K}) where K

Utility function that removes all type-specifiers listed in the curly brackets.

# Examples
```julia-repl
julia> remove_curly(Array{Float64,1})
Array
```
"""
@generated function remove_curly(::Type{K}) where K
    name_without_curly = Meta.parse(string(K)).args[1]
    :( $name_without_curly )
end

"""
    get_curly(::Type{K}) where K

Utility function that returns a tuple with all type-specifiers listed in the
curly brackets.

# Examples
```julia-repl
julia> remove_curly(Array{Float64,1})
(:Float64, 1)
```
"""
@generated function get_curly(::Type{K}) where K
    curly_arguments = Meta.parse(string(K)).args[2:end]
    :( Tuple($curly_arguments) )
end
