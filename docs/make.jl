using SawyerEliassenSolver
using Documenter

DocMeta.setdocmeta!(
    SawyerEliassenSolver, :DocTestSetup, :(using SawyerEliassenSolver); recursive=true
)

pages = [
    "Home" => "index.md",
    "Library" => [
        "Public" => "library/public.md",
        "Internals" => "library/internals.md"
    ],
    "Examples" => [
        "Barotropic NIWs" => "examples/barotropic_niws.md",
        "Variably Stratified Frontal Zone" => "examples/variably_stratified_frontal_zone.md",
        "Baroclinic Vorticity Filament" => "examples/baroclinic_vorticity_filament.md"
    ]
]

makedocs(;
    modules=[SawyerEliassenSolver],
    authors="Jamie Hilditch <hilditch@stanford.edu> and contributors",
    sitename="SawyerEliassenSolver.jl",
    format=Documenter.HTML(;
        canonical="https://Jamie-Hilditch.github.io/SawyerEliassenSolver.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=pages,
    checkdocs=:exports
)

deploydocs(; repo="github.com/Jamie-Hilditch/SawyerEliassenSolver")
