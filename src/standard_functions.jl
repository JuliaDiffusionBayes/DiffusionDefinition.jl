"""
    dimension(::DiffusionProcess{T,DP,DW})

Return dimension of the stochastic process and driving Brownian motion.
"""
dimension(::DiffusionProcess{T,DP,DW}) where {T,DP,DW} = (
    process = DP,
    wiener = DW
)

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

clone(P::T, args...) where T<:DiffusionProcess = T(args...)

function update_params(P::T, new_params...) where T<:DiffusionProcess
    ei = end_point_info(P)
    ei === nothing ? T(new_params...) : T(new_params..., ei...)
end


custom_zero(D::Integer, ::Type{K}) where K <: Array = zeros(eltype(K), D)
custom_zero(D::Integer, ::Type{K}) where K = zero(K)
