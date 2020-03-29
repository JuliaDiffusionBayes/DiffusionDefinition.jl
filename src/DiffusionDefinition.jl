module DiffusionDefinition

    using Random, Trajectories
    using LinearAlgebra, StaticArrays
    using MacroTools

    const ℝ{N} = SVector{N,Float64}

    include("types.jl")
    include("state_space_restrictions.jl")
    include("diffusion_process.jl")
    include(joinpath("..", "examples", "example_list.jl"))
    include("buffer.jl")
    include("sampling.jl")

    export @diffusion_process
    export @load_diffusion

    export UnboundedStateSpace, LowerBoundedStateSpace, UpperBoundedStateSpace,
        BoundedStateSpace

    export AbstractBuffer
    export Wiener

    using Reexport
    @reexport using Trajectories
end # module
