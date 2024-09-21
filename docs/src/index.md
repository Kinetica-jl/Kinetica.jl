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

### Other Kinetica Packages

This process is the same for other Kinetica extension packages, such as KineticaKPM.jl, which can be installed in the same way, e.g.

```julia
using Pkg
Pkg.add("KineticaKPM")
```

### Python Dependencies

Kinetica makes use of Python packages such as [RDKit](https://github.com/rdkit/rdkit) and [Open Babel](https://github.com/openbabel/openbabel) internally for extracting information from molecular geometries. Installation of these packages is handled automatically thanks to [CondaPkg.jl](https://github.com/JuliaPy/CondaPkg.jl), which creates a `conda` environment that is isolated to the current project and the packages within it.

This `conda` environment is composable at runtime, so if you were to have both Kinetica.jl and KineticaKPM.jl in the same Julia project, the Python dependencies of both packages would be automatically assembled into a dedicated `conda` environment that is then used internally by both packages through [PythonCall.jl](https://github.com/JuliaPy/PythonCall.jl). 

### xTB

While not a direct dependency, some parts of Kinetica's CRN exploration routines (which act through the [CDE](https://github.com/HabershonLab/cde) code) require an electronic structure code to perform geometry optimisations and energy calculations. Since only approximate geometries and energies are required within CDE, we recommend using the GFN2-xTB method within the [Extended Tight-Binding (xTB) package](https://github.com/grimme-lab/xtb) by Bannwarth et. al. 

An xTB package is included with Kinetica.jl's Python dependencies and will be installed automatically. While not available on your system's PATH by default, it will be accessible to Kinetica. If you wish to also use Kinetica's installed xTB binary outside of Kinetica, an alias can be created (assuming you are in your Julia project's directory) with:

```bash
alias xtb="$(julia --project -e 'using CondaPkg; print(CondaPkg.which("xtb"))')"
```

### Graphviz

[Graphviz](https://graphviz.org/) is an open source software suite for graph visualisation. Kinetica provides an interface to Graphviz through Catalyst.jl (see [Results Analysis](@ref)) and bundles it as a Python dependency, the same as xTB above.

The Graphviz executables are similarly not added to your system's PATH but are available to Kinetica, and the same aliasing procedure as above can be applied if access to these executables is desired outside of Kinetica.

!!! note "Why not Graphviz_jll?"
    Some readers may note that Graphviz is also distributed as a JLL package through the Julia package manager, and installing it this way may be simpler than treating it as a Python dependency. However, [Graphviz_jll.jl](https://github.com/JuliaBinaryWrappers/Graphviz_jll.jl) is currently compiled without some optional dependencies such as GTS, making it much less useful for graphing large CRNs. We therefore fall back to the Conda package for the time being.

## Citing Kinetica

If you use any of the Kinetica packages in your work, please cite the following:

### Kinetica.jl

> ```Gilkes, J., Storr, M. T., Maurer, R. J., & Habershon, S. (2024). Predicting Long-Time-Scale Kinetics under Variable Experimental Conditions with Kinetica.jl. Journal of Chemical Theory and Computation, 20(12), 5196–5214. https://doi.org/10.1021/acs.jctc.4c00333```

### KineticaKPM.jl

> ```Paper coming soon!```