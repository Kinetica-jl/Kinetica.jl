# Saving & Loading

Kinetica allows for saving and loading [`ODESolveOutput`](@ref) objects (see [Results Analysis](@ref)) in a binary JSON (BSON) format using the [BSON.jl](https://github.com/JuliaIO/BSON.jl) package. This allows for generated CRNs and the parameters and results of kinetic simulations to be serialised and efficiently stored, while also being usable in new Julia sessions and retaining the analysis tools detailed in [Results Analysis](@ref).

## Saving

When performing a CRN exploration, saving is handled automatically by passing a directory to the `savedir` keyword argument of [`explore_network`](@ref), as was demonstrated in [Getting Started](@ref). However, [`ODESolveOutput`](@ref)s can also be saved manually using the [`save_output`](@ref) function:

```julia
# res = solve_network(solvemethod, sd, rd)
save_output(res, "/path/to/saved_output.bson")
```

!!! note "Save Path"
    Unlike the `savedir` argument of [`explore_network`](@ref), the path provided to [`save_output`](@ref) should be to a file ending in `.bson`, rather than to a directory, since only a single CRN is being saved (rather than the multiple checkpoints saved by iterative CRN explorations).

Saving simulation outputs this way destructures them into core Julia objects only - mostly `Vector`s, `Dict`s and `Array`s. This allow them to be loaded back in under any Julia environment, even those where Kinetica or DifferentialEquations.jl (which [`ODESolveOutput`](@ref)s makes heavy use of types from) are not loaded. However, this does necessitate discarding some of the internals of objects such as [`ConditionSet`](@ref)s and DiffEq's `ODESolution`s which are not easily serialised.

!!! warning "Unsaved Fields"
    Among the data that will be lost on saving are:

    * Anything in `ODESolveOutput.sd.cache`, as the contents of this cache are intentionally undefined and therefore difficult to serialise correctly. This will be replaced with an empty `Dict` upon loading.
    * Any Julia functions within variable condition profiles, i.e. the `f` field within `AbstractDirectProfile`s and the `grad` field within `AbstractGradientProfile`s. Functions are not directly serialisable within BSON.
    * The exact ODE solver used within `ODESolveOutput.pars.solver`, as this is always an external type. This will be replaced with a `Symbol` representing the type of solver for reference of what was used within the attached simulation.
    * Any DiffEq `ODESolution`s such as `ODESolveOutput.sol` will lose their underlying ModelingToolkit representations, as only the `sol.u` and `sol.t` fields will be saved. These are reconstructed as `DiffEqArray`s when loaded back in.
    
    Note that no information about the kinetic calculator used within the simulation is saved either. Calculator implementation varies too much to be consistently serialised so it is intentionally left out of [`ODESolveOutput`](@ref). It can be very useful to save some of the information stored within calculators and this is entirely possible with BSON.jl, but this is left to the user.

## Loading

Saved CRN results can be loaded back in to Julia through two methods. If loaded back in using Kinetica's [`load_output`](@ref) function, the serialised BSON gets reconstructed into a new [`ODESolveOutput`](@ref). Taking the CRN we generated and simulated in [Getting Started](@ref) as an example:

```@example saving_loading
using Kinetica

res = load_output("../my_CRN_out/direct_network_final.bson");
nothing # hide
```

!!! note "Suppressing Output"
    When called in an interactive session, [`load_output`](@ref) will produce a lot of output detailing the [`ODESolveOutput`](@ref) object. To suppress this, the load call can be ended with a semicolon, like above.

While every effort is made to maintain compatibility between saved outputs between Kinetica versions, we cannot guarantee that the internal structure of [`ODESolveOutput`](@ref)s will never change. If this happens in a way that cannot be worked around when reconstructing within [`load_output`](@ref), or if you wish to load the data in an environment without Kinetica, the output can still be loaded in as a raw BSON `Dict` tree:

```@example saving_loading
using BSON
out_raw = BSON.load("../my_CRN_out/direct_network_final.bson")
```

This allows users to flexibly access generated CRNs and kinetic simulation results, even when Kinetica is unavailable.