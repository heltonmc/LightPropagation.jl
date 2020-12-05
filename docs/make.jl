#push!(LOAD_PATH,"/home/heltonmc/git-hub/LightPropagation/src")
push!(LOAD_PATH,"/Users/michaelhelton/LightPropagation/src")

using Documenter
using LightPropagation
#=
makedocs(
    sitename = "LightPropagation.jl",
    doctest = true,
    build = "build",
    clean = true,
    source = "src",
    modules= Module[LightPropagation],
    repo = "github.com/heltonmc/LightPropagation.git",
    pages = [
    "Home" => "index.md"
]
)
=#


makedocs(
    sitename = "LightPropagation.jl",
    modules = [LightPropagation],
    format = Documenter.HTML(),
    pages = [ 
                "Home" => "index.md",
                "Forward Models" => Any[
                    "SDA" => "DA_slab_semiinfgeom.md"
                ]

    ]
)

deploydocs(repo = "https://github.com/heltonmc/LightPropagation.jl.git")
