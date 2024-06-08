using CorePore
using Documenter

DocMeta.setdocmeta!(CorePore, :DocTestSetup, :(using CorePore); recursive=true)

makedocs(;
    modules=[CorePore],
    authors="Graham Harper Edwards, Sarah Neuhaus",
    repo="https://github.com/grahamedwards/CorePore.jl/blob/{commit}{path}#{line}",
    sitename="CorePore.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://grahamedwards.github.io/CorePore.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/grahamedwards/CorePore.jl",
    devbranch="main",
)
