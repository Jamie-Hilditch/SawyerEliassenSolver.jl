using SawyerEliassenSolver
using Documenter

DocMeta.setdocmeta!(SawyerEliassenSolver, :DocTestSetup, :(using SawyerEliassenSolver); recursive=true)

makedocs(;
    modules=[SawyerEliassenSolver],
    authors="Jamie Hilditch <hilditch@stanford.edu> and contributors",
    sitename="SawyerEliassenSolver.jl",
    format=Documenter.HTML(;
        canonical="https://Jamie-Hilditch.github.io/SawyerEliassenSolver.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Jamie-Hilditch/SawyerEliassenSolver.jl",
    devbranch="master",
)
