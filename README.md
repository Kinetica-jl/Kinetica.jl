# KineticaCore

## Installation

When installing *KineticaCore* from Julia's Package Registry, a local binary of CDE will need to be compiled. CDE is included as a submodule in *KineticaCore*'s GitHub reposity, and can be automatically compiled in the correct location when building the package. As such, any environment variables that are required for CDE compilation should be set before instantiating the Julia REPL where *KineticaCore* is installed/built.

For example, if CDE is to be compiled with `ifort`, this should be set within the environment:

```bash
$ FC=ifort julia
julia> using Pkg
julia> Pkg.add("KineticaCore")
```

This will pass the Fortran compiler through to *KineticaCore*'s build file `deps/build.jl`, which will compile CDE.

In case of any errors, full details are located in the build log at `deps/build.log`. While it should cover most setups, CDE's `Makefile` may need to be modified for compilation to succeed. After a failed build (which instantiated the git submodule for CDE), the `Makefile` can be accessed from `submodules/cde/Makefile`. After any edits have been made, *KineticaCore* can be rebuilt from the REPL to recompile CDE:

```bash
$ FC=ifort julia
julia> using Pkg
julia> Pkg.build("KineticaCore")
```
