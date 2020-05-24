#===============================================================================

                                Cloning

                This is quite inefficient, but has a nice syntax

===============================================================================#

"""
    clone(P::T, θ::AbstractDict) where T <: DiffusionProcess

Simplified cloning of diffusion law `P`. Substitute relevant parameters with
new values. `θ` must be a dict corresponding to parameters returned after a call
to `var_parameters`.
"""
function clone(P::T, θ::AbstractDict) where T <: DiffusionProcess
    cp = const_parameters(P)
    remove_curly(T)(;cp..., θ...)
end


#===============================================================================

                        in-place parameter setting

            Efficient way of re-parameterizing the diffusion law

===============================================================================#
"""
    set_parameters!(P::DiffusionProcess, θ, entries)

Set parameters of a diffusion law `P` in-place. `entries` should be a collection
of pairs Pair{Int64,Symbol} that list the relevant entries in `θ` for
reparameterization, together with the corresponding parameter names.
"""
function set_parameters!(P::DiffusionProcess, θ, entries)
    for (i,s) in entries
        setfield!(P, s, θ[i])
    end
end

function set_parameters!(P::DiffusionProcess, P°::DiffusionProcess, entries)
    for e in entries
        setfield!(P°, e, getfield(P, e))
    end
end

function set_parameters!(P::T, θ::Dict) where T <: DiffusionProcess
    for k in keys(θ)
        hasfield(T, k) && setfield!(P, k, θ[k])
    end
end

function same_entries(P::DiffusionProcess, P°::DiffusionProcess, entries)
    for e in entries
        getfield(P, e) == getfield(P°, e) || return false
    end
    true
end

function has_any(P::T, entries) where T <: DiffusionProcess
    for (i,s) in entries
        hasfield(T, s) && return true
    end
    false
end
