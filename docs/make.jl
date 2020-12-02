push!(LOAD_PATH,"/home/heltonmc/git-hub/LightPropagation/src")
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
deploydocs(repo = "https://github.com/heltonmc/LightPropagation.jl.git")


makedocs(
    sitename = "LightPropagation.jl",
    format = Documenter.HTML(),
    pages = Any[ "Home" => "index.md"]
)