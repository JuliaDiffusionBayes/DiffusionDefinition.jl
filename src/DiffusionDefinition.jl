module DiffusionDefinition

include("types.jl")
include("state_space_restrictions.jl")
include("diffusion_process.jl")

_DIR = "../examples"
include(joinpath(_DIR, "lorenz_system", "lorenz_system.jl"))
include(joinpath(_DIR, "lotka_volterra", "lotka_volterra.jl"))
include(joinpath(_DIR, "lotka_volterra", "lotka_volterra_aux.jl"))

export @diffusion_process

export UnboundedStateSpace, LowerBoundedStateSpace, UpperBoundedStateSpace,
    BoundedStateSpace

end # module
