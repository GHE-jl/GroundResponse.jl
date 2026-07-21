using GroundResponse
using Documenter

# Make `using GroundResponse` available to every doctest in docstrings and pages.
DocMeta.setdocmeta!(
    GroundResponse,
    :DocTestSetup,
    :(using GroundResponse);
    recursive = true,
)

makedocs(;
    modules = [GroundResponse],
    authors = "Gabriel-Dion <dion.gabriel100@gmail.com>",
    sitename = "GroundResponse.jl",
    format = Documenter.HTML(;
        canonical = "https://GHE-jl.github.io/GroundResponse.jl",
        edit_link = "main",
        assets = String[],
        mathengine = Documenter.KaTeX(),
        sidebar_sitename = false,
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Modeling theory" => [
            "Overview" => "theory/overview.md",
            "Line-source models" => "theory/line_source.md",
            "Moving-source models" => "theory/moving_source.md",
            "Spatial superposition" => "theory/superposition.md",
        ],
        "Borefields" => "borefield.md",
        "API reference" => "api.md",
        "References" => "references.md",
    ],
    # Keep the build strict so broken cross-references or missing docstrings fail CI.
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/GHE-jl/GroundResponse.jl",
    devbranch = "main",
)
