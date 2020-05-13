using Documenter, DiffusionDefinition

makedocs(;
    modules=[DiffusionDefinition],
    format=Documenter.HTML(
        mathengine = Documenter.MathJax(
            Dict(
                :TeX => Dict(
                    :equationNumbers => Dict(
                        :autoNumber => "AMS"
                    ),
                    :Macros => Dict(
                        :dd => "{\\textrm d}",
                        :RR => "\\mathbb{R}",
                        :wt => ["\\widetilde{#1}", 1]
                    ),
                )
            )
        ),
        collapselevel = 1,
    ),
    pages=[
        "Home" => "index.md",
        "Get started" => joinpath("get_started", "package_overview.md"),
        "User manual" => Any[
            "Defining diffusions" => joinpath("manual", "defining_diffusions.md"),
            "Drift and volatility" => joinpath("manual", "drift_and_volatility.md"),
            "Sampling trajectories" => joinpath("manual", "sampling.md"),
            "Computing path-functionals" => joinpath("manual", "functionals_of_paths.md"),
            "Loading diffusions" => joinpath("manual", "load_diff.md"),
            "Utility functions" => joinpath("manual", "convenience_functions.md"),
            "Buffers" => joinpath("manual", "buffers.md"),
            "State space restrictions" => joinpath("manual", "state_space.md"),
            "Reparameterizations" => joinpath("manual", "reparameterization.md"),
        ],
        "How to..." => Any[
            "How to sample Brownian bridges" => joinpath("how_to_guides", "sample_brownian_bridges.md"),
            "(TODO) How to sample Diffusion bridges" => joinpath("how_to_guides", "sample_diffusion_bridges.md"),
            "How to sample complex-valued diffusions" => joinpath("how_to_guides", "sample_complex_valued.md"),
            "(TODO) How to customize my diffusion plots" => joinpath("how_to_guides", "customize_my_diffusion_plots.md"),
            "(TODO) How to compute expectation of gradients of functionals" => joinpath("how_to_guides", "expect_of_gradient.md"),
            "(TODO) How to compute higher order derivatives of functionals" => joinpath("how_to_guides", "higher_order_derivatives.md"),
        ],
        "Tutorials" => Any[
            "(TODO) Sampling Lorenz system" => joinpath("tutorials", "simple_lorenz.md"),
            "(TODO) In-place vs out-of-place methods" => joinpath("tutorials", "in_place_vs_out_of_place.md"),
        ],
        "Predefined Diffusions" => Any[
            "Lotka-Volterra system" => joinpath("predefined_processes", "lotka_volterra.md"),
            "Sine diffusion" => joinpath("predefined_processes", "sine.md"),
            "SIR model" => joinpath("predefined_processes", "sir.md"),
            "Lorenz63 system" => joinpath("predefined_processes", "lorenz63.md"),
            "Lorenz96 system" => joinpath("predefined_processes", "lorenz96.md"),
            "Favetto-Samson model" => joinpath("predefined_processes", "favetto_samson.md"),
            "Prokaryotic autoregulatory gene network" => joinpath("predefined_processes", "prokaryote.md"),
            "FitzHugh-Nagumo model" => joinpath("predefined_processes", "fitzhugh_nagumo.md"),
            "Jansen-Rit model" => joinpath("predefined_processes", "jansen_rit.md"),
        ],
        "Index" => "post_index.md",
    ],
    repo="https://github.com/JuliaDiffusionBayes/DiffusionDefinition.jl/blob/{commit}{path}#L{line}",
    sitename="DiffusionDefinition.jl",
    authors="Sebastiano Grazzi, Frank van der Meulen, Marcin Mider, Moritz Schauer",
    #assets=String[],
)

deploydocs(;
    repo="github.com/JuliaDiffusionBayes/DiffusionDefinition.jl",
)
