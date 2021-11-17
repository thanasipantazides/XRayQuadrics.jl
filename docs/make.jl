using XRayQuadrics
using Documenter

DocMeta.setdocmeta!(XRayQuadrics, :DocTestSetup, :(using XRayQuadrics); recursive=true)

makedocs(;
    modules=[XRayQuadrics],
    authors="Thanasi Pantazides <thanasipantazides@gmail.com> and contributors",
    repo="https://github.com/thanasipantazides/XRayQuadrics.jl/blob/{commit}{path}#{line}",
    sitename="XRayQuadrics.jl",
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