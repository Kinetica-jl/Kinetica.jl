# Kinetica.jl

[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

Kinetica.jl is the core package of the Kinetica organization, used for automated exploration and time integration of large chemical reaction networks (CRNs).

By building on packages within the Julia language's SciML organization (namely [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) and [Catalyst.jl](https://github.com/SciML/Catalyst.jl)), Kinetica enables fast, adaptable CRN exploration, construction and solution under arbitrary variable simulation conditions.

Kinetica provides automated routines for guided CRN exploration, where chemical reactions are only explored if they are predicted to be relevant to the kinetics imposed by the simulation conditions. This requires repeated kinetic modelling of CRNs as they are built, which is facilitated by a discrete approximation to variable rate constant kinetics.

> [!NOTE]
> Kinetica v0.7 updates many dependencies, including [StableHashTraits](https://github.com/beacon-biosignals/StableHashTraits.jl), which it uses for ensuring the uniqueness of reactions in a CRN. CRNs saved in BSON format under previous Kinetica versions may become unusable without re-hashing all their reactions, which can be achieved with the `get_rhash` function. Raw CRNs still in directory tree format (usually imported with `import_network`) are unaffected, as hashing occurs after this stage.

## Documentation

For information on installation, usage and development of Kinetica.jl, see the [documentation](https://kinetica-jl.github.io/Kinetica.jl/stable/).

## Citation

If you use Kinetica.jl in your research, please cite the following paper:

> ```Gilkes, J., Storr, M. T., Maurer, R. J., & Habershon, S. (2024). Predicting Long-Time-Scale Kinetics under Variable Experimental Conditions with Kinetica.jl. Journal of Chemical Theory and Computation, 20(12), 5196–5214. https://doi.org/10.1021/acs.jctc.4c00333```

## License

This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg