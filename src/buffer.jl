#==============================================================================#
#
#          Structs defining buffers used for in-place computations
#
#==============================================================================#

"""
    AbstractBuffer{T} <: AbstractArray{T,1}

Types inheriting from `AbstractBuffer` define buffers used for various in-place
computations. Each subtype MUST:
- have a field `data` which is a one-dimensional vector with data
- provide a constructor of the form: NAME{eltype}(dim1, dim2, ...)
and it CAN:
- have multiple fields that provide views or reshaped views to `data`
"""
abstract type AbstractBuffer{T} <: AbstractArray{T,1} end

const _AB{T} = AbstractBuffer{T} where T # alias

Base.size(var::_AB) = size(var.data)
Base.eltype(var::_AB{T}) where T = T

Base.getindex(var::_AB, i::Int) = var.data[i]
Base.getindex(var::_AB, I::Vararg{Int,N}) where N = var.data[I...]
Base.getindex(var::_AB, ::Colon) = var.data[:]
Base.getindex(var::_AB, kr::AbstractRange) = var.data[kr]

Base.setindex!(var::_AB, v, i::Int) = (var.data[i] = v)
Base.setindex!(var::_AB, v, I::Vararg{Int,N}) where N = (var.data[I...] = v)
Base.setindex!(var::_AB, v, ::Colon) = (var.data[:] .= v)
Base.setindex!(var::_AB, v, kr::AbstractRange) = (var.data[kr] .= v)

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

"""
    dimensions(::Type{K}) where K<:_AB

Utility function for structs inheriting from AbstractBuffer. Returns the
dimensions required for constructing the internal views of the data.
"""
dimensions(::Type{K}) where K<:_AB = get_curly(K)[2].args

Base.similar(d::K) where K <: _AB = begin
    new_data = similar(d.data)
    K(new_data, new_data.x...)
end
# the one below is likely wrong, but it's also not used atm, come back later
Base.similar(d::K, ::Type{T}) where {K <: _AB,T} = begin
    new_data = similar(d.data, T)
    K(new_data, new_data.x...)
end


"""
    struct StandardEulerBuffer{T,D,Tb,Tσ,Tdw} <: AbstractBuffer{T}
        data::Vector{T}
        b::Tb
        σ::Tσ
        dW::Tdw
    end

Standard buffer for the Euler-Maruyama simulations. The data is stored in
`data` and `b`, `σ`, `dW` provide appropriately reshaped views to the
corresponding segments of `data`. `T` is the DataType of each data element
and `D` are dimensions needed to create the views.
"""
struct StandardEulerBuffer{T,Tb,Tσ,Tdw} <: AbstractBuffer{T}
    data::ArrayPartition{T,Tuple{Tb, Tb, Tσ, Tdw}}
    b::Tb
    y::Tb
    σ::Tσ
    dW::Tdw
end

function StandardEulerBuffer(b::Tb, σ::Tσ, dW::Tdw) where {Tb, Tσ, Tdw}
    #TODO fix T for ForwardDiff
    T = eltype(dW)
    @assert eltype(b) == eltype(σ) == T
    y_temp = deepcopy(b)
    StandardEulerBuffer{T, Tb, Tσ, Tdw}(
        ArrayPartition(b, y_temp, σ, dW), b, y_temp, σ, dW
    )
end

function StandardEulerBuffer{K}(P::DiffusionProcess{T,DP}) where {K,T,DP}
    dW = zero(P, Val(:wiener))
    b = zero(K, DP, ismutable(K))
    σ = init_σ_for_buffer(P, Val(diagonaldiff(P)), Val(sparsediff(P)), K)
    StandardEulerBuffer(b, σ, dW)
end

function init_σ_for_buffer(
        P::DiffusionProcess{T,DP,DW},
        ::Val{false},
        ::Val{false},
        ::Type{K}
    ) where {T,DP,DW,K}
    σ = zero(K, (DP, DW), ismutable(K))
end

function init_σ_for_buffer(
        P::DiffusionProcess{T,DP,DW},
        ::Val{true},
        ::Any,
        ::Type{K}
    ) where {T,DP,DW,K}
    σ = Diagonal{eltype(K)}(zero(K, DW, ismutable(K)))
end

function init_σ_for_buffer(
        P::DiffusionProcess{T,DP,DW},
        ::Val{false},
        ::Val{true},
        ::Type{K}
    ) where {T,DP,DW,K}
    error("sparse arrays not implemented yet")
end

"""
    struct LinearDiffBuffer{T,D,Tb,TB,Tσ,Tdw} <: AbstractBuffer{T}
        data::Vector{T}
        b::Tb
        B::TB
        σ::Tσ
        dW::Tdw
    end

A buffer for Euler-Maruyama simulations of linear diffusions. Almost the same
as `StandardEulerBuffer`, but contains additional space for an intermediate
construction of a matrix `B` (and a corrsponding view).
"""
struct LinearDiffBuffer{T,Tb,Tσ,Tdw,TB} <: AbstractBuffer{T}
    data::ArrayPartition{T,Tuple{Tb, Tb, Tσ, Tdw, TB}}
    b::Tb
    y::Tb
    σ::Tσ
    dW::Tdw
    B::TB
end

function LinearDiffBuffer(b::Tb, σ::Tσ, dW::Tdw, B::TB) where {Tb, Tσ, Tdw, TB}
    #TODO fix T for ForwardDiff
    T = eltype(dW)
    @assert eltype(b) == eltype(σ) == eltype(TB) == T
    y_temp = deepcopy(b)
    LinearDiffBuffer{T, Tb, Tσ, Tdw, TB}(
        ArrayPartition(b, y_temp, σ, dW, B), b, y_temp, σ, dW, B
    )
end

function LinearDiffBuffer{K}(P::LinearDiffusion{T,DP,DW}) where {K,T,DP,DW}
    dW = zero(K, DW, ismutable(K))
    b = zero(K, DP, ismutable(K))
    σ = init_σ_for_buffer(P, Val(diagonaldiff(P)), Val(sparsediff(P)), K)
    B = zero(K, (DP, DP), ismutable(K)) # change to initialize similarly to σ
    LinearDiffBuffer(b, σ, dW, B)
end

get_B(seb::LinearDiffBuffer) = seb.data.x[5]
