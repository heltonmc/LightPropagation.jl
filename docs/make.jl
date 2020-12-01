push!(LOAD_PATH,"/home/heltonmc/git-hub/LightPropagation/src")
using Documenter
using LightPropagation

makedocs(
    sitename = "LightPropagation.jl",
    doctest = false,
    modules=[LightPropagation],
    repo = "github.com/heltonmc/LightPropagation.git",
    pages = [
    "Home" => "index.md",
    "DT"   => "DAslab_semiinfgeom.md"
]
)

deploydocs(repo = "https://github.com/heltonmc/LightPropagation.jl.git")