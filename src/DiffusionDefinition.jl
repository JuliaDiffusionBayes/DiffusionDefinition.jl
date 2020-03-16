module DiffusionDefinition

include("types.jl")
include("state_space_restrictions.jl")
include("diffusion_process.jl")

export @diffusion_process

export UnboundedStateSpace, LowerBoundedStateSpace, UpperBoundedStateSpace,
    BoundedStateSpace

end # module
