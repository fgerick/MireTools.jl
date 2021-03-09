using MireTools
using Documenter

DocMeta.setdocmeta!(MireTools, :DocTestSetup, :(using MireTools); recursive=true)

makedocs(;
    modules=[MireTools],
    authors="fgerick <felixgerick@gmail.com> and contributors",
    repo="https://github.com/fgerick/MireTools.jl/blob/{commit}{path}#{line}",
    sitename="MireTools.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fgerick.github.io/MireTools.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fgerick/MireTools.jl",
)
