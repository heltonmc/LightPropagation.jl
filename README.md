# LightPropagation

## Purpose

LightPropagation provides a set of tools to model and analyze the propagation of light in turbid media written in the [Julia programming language](https://julialang.org/).
The current library supports simulating the time-resolved reflectance in the forward direction under the diffusion approximation for 3 different geometries:
- __Semi-infinite__ 
- __Slabs__
- __Parallelepipeds__

## Usage

The easiest way to install julia is by [downloading julia](https://julialang.org/downloads/) from the offical site and following the [platform specific installations](https://julialang.org/downloads/platform/). 

The current package is not registered with the general Julia repository. 
The easiest way to use is to clone the repository by typing `git clone https://github.com/heltonmc/LightPropagation.git` in your terminal.

Launch julia and in the REPL cd into the LightPropagation folder you just cloned. You can check your current working directory using `pwd()` in the julia REPL.
```julia
cd("...\LightPropagation")
] activate .
using LightPropagation
```
Once your current directory is the LightPropagation repository you can activate it as shown above. In julia you can activate the package manager by typing `]` into the julia repl and then typing `activate .`. You can also use shell commands by typing `;` into the repl. To go back to the julia environment simplay backspace.
To load packages in julia type `using LightPropagation` in the julia repl.

## Tutorials

## Documentation

## LightPropagation community

Currently none.

## Contributing to LightPropagation

LightPropgation is a collaborative project open to contributions and discussion.

## How to cite LightPropagation

In order to give credit to the `LightPropagation` contributors, we ask you to cite the reference below in any publication in which you have made use of `LightPropagation` packages:

## Contact

Please contact the project adminstrators, [Michael Helton](mailto:heltonmc@umich.edu), for further questions on usage and contributing. 

