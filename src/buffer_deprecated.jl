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
