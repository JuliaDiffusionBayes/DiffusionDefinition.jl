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

const _AB{T} = AbstractBuffer{T} where T # an alias

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



@generated function Base.similar(::K) where K<:_AB{T} where T
    pure_name = Meta.parse(string((remove_curly(K))))
    constructor = Expr(
        :call,
        Expr(
            :curly,
            :($pure_name),
            T,
        )
    )
    append!(constructor.args, [dimensions(K)...])
    :(
        $constructor
    )
end

@generated function Base.similar(::K, ::Type{T}) where {K<:_AB,T}
    pure_name = Meta.parse(string((remove_curly(K))))
    constructor = Expr(
        :call,
        Expr(
            :curly,
            :($pure_name),
            T,
        )
    )
    append!(constructor.args, [dimensions(K)...])
    :(
        $constructor
    )
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
struct StandardEulerBuffer{T,D,Tb,Tσ,Tdw} <: AbstractBuffer{T}
    data::Vector{T}
    b::Tb
    σ::Tσ
    dW::Tdw

    function StandardEulerBuffer{T}(dim_process::Number, dim_wiener::Number) where T
        D, m = dim_process, dim_wiener
        data = zeros(T, D + D*m + m)

        b = view(data, 1:D)
        σ = reshape( view(data, (D+1):(D+D*m)), (D,m) )
        dW = view(data, (D+D*m+1):(D+D*m+m))
        new{T,(D,m),typeof(b),typeof(σ),typeof(dW)}(data, b, σ, dW)
    end
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
struct LinearDiffBuffer{T,D,Tb,TB,Tσ,Tdw} <: AbstractBuffer{T}
    data::Vector{T}
    b::Tb
    B::TB
    σ::Tσ
    dW::Tdw

    function LinearDiffBuffer{T}(dim_process::Number, dim_wiener::Number) where T
        D, m = dim_process, dim_wiener
        data = zeros(T, D + D*D + D*m + m)

        b = view(data, 1:D)
        B = reshape( view(data, (D+1):(D+D*D)), (D,D) )
        σ = reshape( view(data, (D+D*D+1):(D+D*D+D*m)), (D,m) )
        dW = view(data, (D+D*D+D*m+1):(D+D*D+D*m+m))
        new{T,(D,m),typeof(b),typeof(B),typeof(σ),typeof(dW)}(data, b, B, σ, dW)
    end
end
