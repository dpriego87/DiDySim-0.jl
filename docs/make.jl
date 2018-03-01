using Documenter, DiDySim

makedocs(;
    modules=[DiDySim],
    format=:html,
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/dpriego87/DiDySim.jl/blob/{commit}{path}#L{line}",
    sitename="DiDySim.jl",
    authors="Daniel Priego-Espinosa",
    assets=[],
)

deploydocs(;
    repo="github.com/dpriego87/DiDySim.jl",
    target="build",
    julia="0.6",
    deps=nothing,
    make=nothing,
)
