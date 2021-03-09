using KerrQuasinormalModes
using Documenter

DocMeta.setdocmeta!(KerrQuasinormalModes, :DocTestSetup, :(using KerrQuasinormalModes); recursive=true)

makedocs(;
    modules=[KerrQuasinormalModes],
    authors="Acme Corp",
    repo="https://github.com/Potatoasad/KerrQuasinormalModes.jl/blob/{commit}{path}#{line}",
    sitename="KerrQuasinormalModes.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Potatoasad.github.io/KerrQuasinormalModes.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Potatoasad/KerrQuasinormalModes.jl",
)
