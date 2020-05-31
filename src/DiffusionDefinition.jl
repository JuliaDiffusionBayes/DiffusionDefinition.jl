module DiffusionDefinition

    using Random, Trajectories
    using LinearAlgebra, StaticArrays, SparseArrays
    using MacroTools
    using RecipesBase
    using RecursiveArrayTools
    using ForwardDiff
    import ForwardDiff: Dual, Tag

    const ℝ{N} = SVector{N,Float64}

    include("types.jl")
    include("standard_functions.jl")
    include("state_space_restrictions.jl")
    include(joinpath("..", "examples", "example_list.jl"))
    include("diffusion_process.jl")
    include("buffer.jl")
    include("trajectories_extensions.jl")
    include("sampling.jl")
    include("reparameterizations.jl")
    include("conjugate_updates.jl")

    export @diffusion_process
    export @load_diffusion
    export @load_variable_diffusion
    export @conjugate_gaussian

    export UnboundedStateSpace, LowerBoundedStateSpace, UpperBoundedStateSpace,
        BoundedStateSpace

    export AbstractBuffer
    export Wiener, wiener

    using Reexport
    @reexport using Trajectories
end # module
