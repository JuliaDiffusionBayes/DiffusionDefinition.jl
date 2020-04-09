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
        )
    ),
    pages=[
        "Home" => "index.md",
        "Overview" => Any[
            "Defining diffusions" => joinpath("overview", "defining_diffusions.md"),
            "Loading diffusions" => joinpath("overview", "load_diff.md"),
            "Sampling trajectories" => joinpath("overview", "sampling.md"),
        ],
        "Examples" => Any[
            "Lotka-Volterra system" => joinpath("examples", "lotka_volterra.md"),
            "Sine diffusion" => joinpath("examples", "sine.md"),
            "SIR model" => joinpath("examples", "sir.md"),
            "Lorenz63 system" => joinpath("examples", "lorenz63.md"),
            "Lorenz96 system" => joinpath("examples", "lorenz96.md"),
            "Favetto-Samson model" => joinpath("examples", "favetto_samson.md"),
            "Prokaryotic autoregulatory gene network" => joinpath("examples", "prokaryote.md"),
            "FitzHugh-Nagumo model" => joinpath("examples", "fitzhugh_nagumo.md"),
            "Jansen-Rit model" => joinpath("examples", "jansen_rit.md"),
        ],
        "Index" => "post_index.md",
    ],
    repo="https://github.com/JuliaDiffusionBayes/DiffusionDefinition.jl/blob/{commit}{path}#L{line}",
    sitename="DiffusionDefinition.jl",
    authors="Sebastiano Grazzi, Frank van der Meulen, Marcin Mider, Moritz Schauer",
    assets=String[],
)

deploydocs(;
    repo="github.com/JuliaDiffusionBayes/DiffusionDefinition.jl",
)
