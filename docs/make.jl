using Documenter
using LightPropagation

makedocs(modules = [LightPropagation], sitename = "LightPropagation.jl")
deploydocs(repo = "https://github.com/heltonmc/LightPropagation.jl.git")