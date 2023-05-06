using CEFOracle
using Documenter

DocMeta.setdocmeta!(CEFOracle, :DocTestSetup, :(using CEFOracle); recursive=true)

makedocs(;
    modules=[CEFOracle],
    authors="Abe Hernández",
    repo="https://github.com/abehersan/CEFOracle.jl/blob/{commit}{path}#{line}",
    sitename="CEFOracle.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://abehersan.github.io/CEFOracle.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/abehersan/CEFOracle.jl",
    devbranch="main",
)
