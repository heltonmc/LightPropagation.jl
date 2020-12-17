#push!(LOAD_PATH,"/home/heltonmc/git-hub/LightPropagation/src")
#push!(LOAD_PATH,"/Users/michaelhelton/LightPropagation/src")

using LightPropagation

import Pkg; Pkg.add("Documenter")
using Documenter
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
    modules = [LightPropagation],
    sitename = "LightPropagation.jl",
    format = Documenter.HTML(),
    pages = [ 
                "Home" => "index.md",
                "Forward Models" => Any[
                    "SDA" => "DA_slab_semiinfgeom.md"
                ]

    ]
)

deploydocs(
    repo = "heltonmc.github.io/LightPropagation",
    devbranch = "main"
)
