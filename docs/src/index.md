# Kinetica.jl Documentation

[Kinetica.jl](https://github.com/Kinetica-jl/Kinetica.jl) is a Julia package for performing automated exploration of chemical reaction networks (CRNs) and integrating these networks in time. In particular, it features:

### Arbitrary simulation conditions

Kinetica.jl is built around giving users complete freedom over kinetic simulations. As such any combination of customisable simulation conditions, static or variable, can be utilised by binding symbolic variable names to flexible parametric condition profiles.

### Kinetics-driven CRN exploration

Chemical space exploration can be performed with a fully random approach, sampling every reaction within a defined number of intermediates from a starting system. This can be very difficult to sample completely and is often inefficient. We provide a focused kinetics-driven approach that selectively explores reaction space in places relevant to the given simulation conditions.

### Flexible kinetic simulation

By leveraging packages from Julia's [SciML](https://sciml.ai/) organization (including [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) and [Catalyst.jl](https://github.com/SciML/Catalyst.jl)), users can perform difficult long-timescale integrations of generated CRNs under challenging variable experimental conditions. 

We supplement this with a discrete approximation to variable rate constant simulations that greatly improves overall solution efficiency and allows for previously inaccessible levels of theory to be incorporated into variable kinetic calculations.

### Modular kinetic calculators

Extending user control over kinetic simulations, Kinetica.jl makes use of a modular calculator interface for rate constant calculations. This allows for a wide variety of techniques to be utilised within kinetic simulations, ranging from expensive DFT-based approaches to fast ML-based approximations.

We currently provide the [KineticaKPM.jl](https://github.com/Kinetica-jl/KineticaKPM.jl) package for calculating rate constants from ML-predicted activation energies, and aim to release further calculator packages in the future. However, the calculator interface allows for simple user definition of new methods too.

## Installation

Kinetica.jl can be installed through the Julia package manager by adding the KinetcaRegistry package registry:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url="https://github.com/Kinetica-jl/KineticaRegistry"))
Pkg.add("Kinetica")
```

This will fetch the latest version of Kinetica.jl, as well as all of its dependencies.

While Kinetica is a Julia package, it has numerous Python-based dependencies that are called through [PyCall.jl](https://github.com/JuliaPy/PyCall.jl), such as Open Babel and RDKit. These have complex multi-language dependencies of their own, so they are most easily installed through a package manager such as [Conda.jl](https://github.com/JuliaPy/Conda.jl). To ensure PyCall correctly uses Conda as its Python environment, the `PYTHON` environment variable should be blank (see below).

By default, all Kinetica packages will install their Python dependencies within their build phase. This occurs automatically when a package is first installed, but can be triggered manually with the following:

```julia
ENV["PYTHON"] = ""
Pkg.build("Kinetica")
```

Manually triggering Kinetica's build phase also does so for PyCall.jl and Conda.jl, so running this will correctly link Conda's Python environment to PyCall, and also install all of Kinetica's Python dependencies.

If Conda-based dependency installation is NOT required (e.g. a rebuild is being performed and the correct Conda environment is already present), this behaviour can be disabled by setting the `KINETICA_BUILD_IGNORE_CONDA` environment variable to `TRUE`, i.e.

```bash
KINETICA_BUILD_IGNORE_CONDA=TRUE julia
julia> using Pkg
julia> Pkg.build("Kinetica")
```

!!! warning "This may take a while..."
    Installing Kinetica's Python dependencies has been known to take some time, so if it looks like the build phase is hanging here, just let it go for a bit longer. If in doubt, you can check the build log to see what Conda is doing (the build phase should output its location).

### Other Kinetica Packages

This process is the same for other Kinetica extension packages, such as KineticaKPM.jl, which can be installed in the same way, e.g.

```julia
using Pkg
Pkg.add("KineticaKPM")
```

If the package being installed has Python dependencies, these should also be installed to the Conda environment in its build phase. Similarly, if this behaviour is not desired, it is also disabled by the `KINETICA_BUILD_IGNORE_CONDA` environment variable, as above.

### xTB

While not a direct dependency, some parts of Kinetica's CRN exploration routines (which act through the [CDE](https://github.com/HabershonLab/cde) code) require an electronic structure code to perform geometry optimisations and energy calculations. Since only approximate geometries and energies are required within CDE, we recommend using the GFN2-xTB method within the [Extended Tight-Binding (xTB) package](https://github.com/grimme-lab/xtb) by Bannwarth et. al. Instructions for installation are available within their [documentation](https://xtb-docs.readthedocs.io/en/latest/setup.html).

!!! note "Note for tutorial users"
    If you plan to follow the tutorials in this documentation, the example CDE inputs assume that you have xTB installed and available system-wide (in your `PATH`). While this can be easily modified to use other electronic structure methods, this is out of the scope of the tutorials and will not be covered. Proceed without xTB installed at your own peril!

## Citing Kinetica

If you use any of the Kinetica packages in your work, please cite the following:

### Kinetica.jl

```bibtex
Paper coming soon!
```

### KineticaKPM.jl

```bibtex
Paper coming soon!
```