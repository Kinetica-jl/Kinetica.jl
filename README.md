# Kinetica.jl

Kinetica.jl is the core package of the Kinetica organization, used for automated exploration and time integration of large chemical reaction networks (CRNs).

By building on packages within the Julia language's SciML organization (namely [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) and [Catalyst.jl](https://github.com/SciML/Catalyst.jl)), Kinetica enables fast, adaptable CRN exploration, construction and solution under arbitrary variable simulation conditions.

Kinetica provides automated routines for guided CRN exploration, where chemical reactions are only explored if they are predicted to be relevant to the kinetics imposed by the simulation conditions. This requires repeated kinetic modelling of CRNs as they are built, which is facilitated by a discrete approximation to variable rate constant kinetics.

## Documentation

For information on installation, usage and development of Kinetica.jl, see the documentation (in progress).

## Citation

If you use Kinetica.jl in your research, please cite the following paper:

```
Citation on the way soon!
```