#==============================================================================#
#
#    Structs responsible for determining the restrictions imposed on the
#    state space of a diffusion process
#
#==============================================================================#

"""
    UnboundedStateSpace <: DiffusionStateSpace

No restrictions imposed on the state-space of the process (i.e. ℝᵈ)
"""
struct UnboundedStateSpace <: DiffusionStateSpace end

"""
    LowerBoundedStateSpace{T,S,N} <: DiffusionStateSpace

Lower bounds imposed on the state-space of a diffusion process. `T` is used
to list the indices that have lower-bound restrictions, `S` indicates the values
of the lower-bounds, `N` is the total number of coordinates with lower-bound
restrictions
"""
struct LowerBoundedStateSpace{T,S,N} <: DiffusionStateSpace

    function LowerBoundedStateSpace(
        coords::NTuple{N,Integer},
        bounds::NTuple{N,Number},
        ) where N
        new{coords, bounds, N}()
    end

    function LowerBoundedStateSpace(
        coords,
        bounds,
        )
        N = length(bounds)
        @assert length(coords) == N
        @assert all( map(c->(typeof(c)<:Integer), coords) )
        @assert all( map(b->(typeof(b)<:Number), bounds) )
        new{Tuple(coords), Tuple(bounds), N}()
    end

end

bound_info(::LowerBoundedStateSpace{T,S,N}) where {T,S,N} = T,S,N

"""
    UpperBoundedStateSpace{T,S,N} <: DiffusionStateSpace

Upper bounds imposed on the state-space of a diffusion process. `T` is used
to list the indices that have upper-bound restrictions, `S` indicates the values
of the upper-bounds, `N` is the total number of coordinates with upper-bound
restrictions
"""
struct UpperBoundedStateSpace{T,S,N} <: DiffusionStateSpace
    function UpperBoundedStateSpace(coords, bounds)
        T,S,N = bound_info(LowerBoundedStateSpace(coords, bounds))
        new{T,S,N}()
    end
end

"""
    BoundedStateSpace{L,U} <: DiffusionStateSpace

Upper and lower bounds imposed on the state-space of a diffusion process.
`L` corresponds to lower bounds, `U` corresponds to upper bounds.
"""
struct BoundedStateSpace{L,U} <: DiffusionStateSpace
    function BoundedStateSpace(
        (coords_lower, bounds_lower),
        (coords_upper, bounds_upper)
        )
        L = LowerBoundedStateSpace(coords_lower, bounds_lower)
        U = UpperBoundedStateSpace(coords_upper, bounds_upper)
        new{L,U}()
    end

    function BoundedStateSpace(
        L::LowerBoundedStateSpace,
        U::UpperBoundedStateSpace
        )
        new{L,U}()
    end
end

"""
    bound_satisfied(::UnboundedStateSpace, x)

No restrictions, bounds satisfied by default
"""
@inline bound_satisfied(::UnboundedStateSpace, x) = true

"""
    bound_satisfied(::LowerBoundedStateSpace{T,S,N}, x) where {T,S,N}

Checks if all coordinates adhere to lower bound restrictions
"""
@generated function bound_satisfied(
    ::LowerBoundedStateSpace{T,S,N},
    x
    ) where {T,S,N}
    ex = :(true)
    for i = 1:N
        ex = :(x[T[$i]] > S[$i] ? $ex : false)
    end
    return ex
end

"""
    bound_satisfied(::UpperBoundedStateSpace{T,S,N}, x) where {T,S,N}

Checks if all coordinates adhere to upper bound restrictions
"""
@generated function bound_satisfied(
    ::UpperBoundedStateSpace{T,S,N},
    x
    ) where {T,S,N}
    ex = :(true)
    for i = 1:N
        ex = :(x[T[$i]] < S[$i] ? $ex : false)
    end
    return ex
end

"""
    bound_satisfied(::BoundedStateSpace{L,U}, x) where {L,U}

Checks if all coordinates adhere to lower and upper bound restrictions
"""
function bound_satisfied(::BoundedStateSpace{L,U}, x) where {L,U}
    bound_satisfied(L, x) && bound_satisfied(U, x)
end

function bound_info(::DiffusionProcess{T,DP,DW,SS}) where {T,DP,DW,SS}
    bound_info(SS)
end

function bound_satisfied(::DiffusionProcess{T,DP,DW,SS}, x) where {T,DP,DW,SS}
    bound_satisfied(SS, x)
end
