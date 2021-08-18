# Installation Instructions

Install Julia by [downloading](https://julialang.org/downloads/) the latest version from the offical site and following the [platform specific installations](https://julialang.org/downloads/platform/).

LightPropagation is not yet a registered Julia package in the official entry, though we still recommend installing it with the built-in Julia package manager. It automatically installs a currently stable and tagged release. From the Julia REPL (accessed by pressing `]` in the Julia REPL command prompt), you can add the package:

```julia
julia>]
(@v1.6) pkg> add "https://github.com/heltonmc/LightPropagation.jl.git"
```

LightPropagation.jl is only tested on Julia v1.5 and later. The package is under rapid development and breaking changes to the user API do occur, so update with care. Please open an issue if something breaks or doesn't seem right!

After adding the package, we can load it with `using` which will load the module and make its exported names availabe for direct use.
```julia
julia> using LightPropagation
```

Throughout the rest of the documenation and examples we will assume that the package has been added and have already typed `using LightPropagation` in the Julia REPL.

### Updating to latest version

If desired, you can also check the version of LightPropagation.jl that you have installed with the `status` command in the package manager...

```julia
julia> ]
(@v1.6) pkg> status LightPropagation
      Status `~/.julia/environments/v1.6/Project.toml`
  [bd080553] LightPropagation v0.4.3 `https://github.com/heltonmc/LightPropagation.jl.git#main`
```

LightPropagation.jl can be updated to the latest tagged release from the package manager by typing...

```julia
(@v1.6) pkg> update LightPropagation
```