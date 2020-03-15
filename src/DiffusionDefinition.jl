module DiffusionDefinition

include("types.jl")
include("state_space_restrictions")
include("diffusion_process.jl")

export @diffusion_process

end # module
