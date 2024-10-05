using SawyerEliassenSolver
using Documenter
using DocumenterCitations
# using DocumenterInterLinks
using Literate

#metadata
DocMeta.setdocmeta!(
    SawyerEliassenSolver, :DocTestSetup, :(using SawyerEliassenSolver); recursive=true
)

##  pages

# theory
sawyer_eliassen_pages = [
    "Theory" => "Sawyer-Eliassen equation/theory.md",
    "Numerics" => "Sawyer-Eliassen equation/numerics.md",
]

# examples
const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR = joinpath(@__DIR__, "src/literated")

example_scripts = [
    "Near-inertial waves" => "near-inertial_waves.jl",
    "Barotropic zeta refraction" => "barotropic_zeta_refraction.jl",
]

example_pages = map(example_scripts) do example
    name = first(example)
    script = last(example)
    example_filepath = joinpath(EXAMPLES_DIR, script)
    Literate.markdown(
        example_filepath, OUTPUT_DIR; flavor=Literate.DocumenterFlavor(), execute=true
    )
    page_name = replace(script, ".jl" => ".md")
    name => joinpath("literated", page_name)
end

# library
library_pages = [
    "SawyerEliassenSolver" => "library/SawyerEliassenSolver.md",
    "Domains" => "library/domains.md",
    "Variables" => "library/variables.md",
    "Forcing" => "library/forcing.md",
    "Problems" => "library/problems.md",
    "Timesteppers" => "library/timesteppers.md",
    "Output" => "library/output.md",
]

pages = [
    "Home" => "index.md",
    "The Sawyer-Eliassen Equation" => sawyer_eliassen_pages,
    "Examples" => example_pages,
    "Library" => library_pages,
    "References" => "references.md",
]

# format
format = Documenter.HTML(;
    prettyurls=get(ENV, "CI", nothing) == "true",
    canonical="https://Jamie-Hilditch.github.io/SawyerEliassenSolver",
    mathengine=MathJax3(),
    assets=String[],
    collapselevel=1,
    size_threshold=500 * 2^10, # 500 KiB
    size_threshold_warn=200 * 2^10, # 200 KiB
)

# references
bib_filepath = joinpath(dirname(@__FILE__), "references.bib")
bib = CitationBibliography(bib_filepath; style=:authoryear)

# interlinks
# links = InterLinks(
#     "HDF5" => "https://juliaio.github.io/HDF5.jl/stable/"
# )

# where the magic happens
makedocs(;
    modules=[SawyerEliassenSolver],
    authors="Jamie Hilditch <hilditch@stanford.edu> and contributors",
    sitename="SawyerEliassenSolver",
    format=format,
    pages=pages,
    plugins=[bib],
    checkdocs=:exports,
)

# clean up h5 files
@info "Cleaning up temporary .h5 created by examples..."
readdir(OUTPUT_DIR; join=true) |> filter(x -> endswith(x, ".h5")) .|> rm

deploydocs(; repo="github.com/Jamie-Hilditch/SawyerEliassenSolver.git")
