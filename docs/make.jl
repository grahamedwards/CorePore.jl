using PorewaterDiffusion
using Documenter

DocMeta.setdocmeta!(PorewaterDiffusion, :DocTestSetup, :(using PorewaterDiffusion); recursive=true)

makedocs(;
    modules=[PorewaterDiffusion],
    authors="Graham Harper Edwards, Sarah Neuhaus",
    repo="https://github.com/grahamedwards/PorewaterDiffusion.jl/blob/{commit}{path}#{line}",
    sitename="PorewaterDiffusion.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://grahamedwards.github.io/PorewaterDiffusion.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/grahamedwards/PorewaterDiffusion.jl",
    devbranch="main",
)
