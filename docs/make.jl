#push!(LOAD_PATH,"/home/heltonmc/git-hub/LightPropagation/src")
#push!(LOAD_PATH,"/Users/michaelhelton/LightPropagation/src")

using LightPropagation

import Pkg; Pkg.add("Documenter")
using Documenter

makedocs(
    modules = [LightPropagation],
    sitename = "LightPropagation.jl",
    format = Documenter.HTML(),
    pages = [ 
                "Home" => "index.md",
                "API" => "API.md",
                "Getting Started" => "getting-started.md",
                "Tutorials" => Any[
                    "When is semi-infinite semi-infinite?" => "whensemiinf.md"
                ],
                "Forward Models" => Any[
                    "Diffusion Approximation" => Any[
                    "Slab" => "DA_slab_semiinfgeom.md",
                    "N-layer Cylinder" => "Nlayer_cyl.md"
                     ],
                ],
    ]
)

deploydocs(
    repo = "github.com/heltonmc/LightPropagation.jl.git",
    devbranch = "main"
)
