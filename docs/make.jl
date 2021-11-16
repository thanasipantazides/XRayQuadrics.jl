using XRayTrace
using Documenter

DocMeta.setdocmeta!(XRayTrace, :DocTestSetup, :(using XRayTrace); recursive=true)

makedocs(;
    modules=[XRayTrace],
    authors="Thanasi Pantazides <thanasipantazides@gmail.com> and contributors",
    repo="https://github.com/thanasipantazides/XRayTrace.jl/blob/{commit}{path}#{line}",
    sitename="XRayTrace.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => [
            "Introduction" => "introduction.md",
            "Quadric Surfaces" => "quadrics.md",
            "Intersections and Collisions" => "collisions.md"
        ]
    ],
)