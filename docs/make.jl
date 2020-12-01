push!(LOAD_PATH,"/home/heltonmc/git-hub/LightPropagation/src")
using Documenter
using LightPropagation

makedocs(
    format = Documenter.HTML(),
    sitename = "LightPropagation.jl",
    doctest = false
)


deploydocs(repo = "https://github.com/heltonmc/LightPropagation.jl.git")