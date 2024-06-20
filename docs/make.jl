using SawyerEliassenSolver
using Documenter
using DocumenterCitations


#metadata
DocMeta.setdocmeta!(
    SawyerEliassenSolver, :DocTestSetup, :(using SawyerEliassenSolver); recursive=true
)

# pages
sawyer_eliassen_pages = [
    "Theory" => "Sawyer-Eliassen equation/theory.md",
    "Numerics" => "Sawyer-Eliassen equation/numerics.md"
]

example_pages = [
    "Barotropic NIWs" => "examples/barotropic_niws.md",
    "Variably Stratified Frontal Zone" => "examples/variably_stratified_frontal_zone.md",
    "Baroclinic Vorticity Filament" => "examples/baroclinic_vorticity_filament.md"
]

library_pages = [
    "Public" => "library/public.md",
    "Internals" => "library/internals.md"
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
    prettyurls = get(ENV, "CI", nothing) == "true",
    canonical="https://Jamie-Hilditch.github.io/SawyerEliassenSolver",
    mathengine = MathJax3(),
    assets=String[],
)

# references
bib_filepath = joinpath(dirname(@__FILE__), "references.bib")
bib = CitationBibliography(bib_filepath, style=:authoryear)

# where the magic happens
makedocs(;
    modules=[SawyerEliassenSolver],
    authors="Jamie Hilditch <hilditch@stanford.edu> and contributors",
    sitename="SawyerEliassenSolver",
    format=format,
    pages=pages,
    plugins=[bib],
    checkdocs=:exports
)

deploydocs(; repo="github.com/Jamie-Hilditch/SawyerEliassenSolver.git")
