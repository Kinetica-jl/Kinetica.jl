# Results Analysis

```@setup results_analysis
mkpath("../assets/tutorials/results_analysis")
using Random
Random.seed!(12345)
```

Upon completion of a kinetic simulation run by [`solve_network`](@ref) (or through the end of a CRN exploration using [`explore_network`](@ref), which also calls this function), Kinetica returns an [`ODESolveOutput`](@ref) object. This binds together the CRN that this simulation was run on (after modifications such as those described in [Filtering CRNs](@ref) and [Removing Low-Rate Reactions](@ref implementation_low_rate)), the parameters the simulation was run with, and the simulation's results. This data can be used for analysing the results of a kinetic simulation.

We will demonstrate this analysis using the results of the CRN we generated and simulated in [Getting Started](@ref):

```@example results_analysis
using Kinetica
res = load_output("../my_CRN_out/direct_network_final.bson")
nothing # hide
```

## Analysing the CRN

The output CRN is stored across two data structures: a [`SpeciesData`](@ref) for holding information about the species within the CRN, and a [`RxData`](@ref) for holding information about the reactions between those species. These can be accessed through `ODESolveOutput.sd` and `ODESolveOutput.rd` respectively. For more information on these objects, see [CRN Representation](@ref crn_representation_page).

```@example results_analysis
sd, rd = res.sd, res.rd
nothing # hide
```

By default, Kinetica stores information about these species and reactions in a compact form to maintain efficiency when working with large CRNs. For example, reactions are stored as `Vector`s of integer species IDs. To inspect individual reactions, we can use the [`print_rxn`](@ref) function:

```@example results_analysis
# Print the first reaction in the CRN.
print_rxn(sd, rd, 1)
nothing # hide
```

!!! note "Using format_rxn"
    [`print_rxn`](@ref) is just a wrapper to print the string returned by [`format_rxn`](@ref), which uses the reactants/products and stoichiometries in `rd` and the species SMILES in `sd` to generate a human-readable string describing a reaction. If you need to print more information alongside a reaction, consider something like `println("$(format_rxn(sd, rd, rid)): $(extra_data[rid])")`.

    Both [`print_rxn`](@ref) and [`format_rxn`](@ref) accept a keyword argument `display_level`, which also tags the requested reaction with the exploration level (see [Iterative CRN Exploration](@ref)) in which it was discovered.

The total number of reactions in a CRN is stored in `RxData.nr`, so printing all the reactions in a CRN is as simple as this loop:

```@example results_analysis
for i in 1:rd.nr
    print_rxn(sd, rd, i)
end
nothing # hide
```

Other fields of [`RxData`](@ref) can also be useful for analysis, such as `RxData.mapped_rxns`. This field contains the atom-mapped reaction SMILES of all reactions in a CRN, and can be used for constructing multi-species reactant and product geometries with consistent atom indices (vital for any reaction path techniques):

```@example results_analysis
println("Atom-mapped reaction SMILES: $(rd.mapped_rxns[1])\n")

# Split reaction SMILES into individual SMILES for reactants/products.
am_reacs, am_prods = String.(split(rd.mapped_rxns[1], ">>"))

# Generate non-overlapping systems of molecules for reactants and products.
reac_species = [sd.toStr[sid] for sid in rd.id_reacs[1]]
reacsys = system_from_smiles(reac_species; dmin=3.0)
prod_species = [sd.toStr[sid] for sid in rd.id_prods[1]]
prodsys = system_from_smiles(prod_species; dmin=3.0)

# Remap atom indices using atom-mapped SMILES and RDKit substructure matching.
mapped_reacs = atom_map_frame(am_reacs, reacsys)
mapped_prods = atom_map_frame(am_prods, prodsys)

# Add extra info to frames, convert to ExtXYZ and print.
mapped_reacs["info"]["type"] = "reactant"
mapped_reacs["info"]["SMILES"] = join(sort(reac_species), ".")
mapped_prods["info"]["type"] = "product"
mapped_prods["info"]["SMILES"] = join(sort(prod_species), ".")
println("Reaction XYZ:")
print(frame_to_xyz(mapped_reacs))
print(frame_to_xyz(mapped_prods))
```

As a result of the atom mapping, the hydrogen atoms which dissociate from the reactant ethane are now consistently indexed between reactant and product geometries.

### Species-Reaction Graphs

Kinetica hooks in to the same graph plotting code as [Catalyst.jl](https://github.com/SciML/Catalyst.jl), which uses [Graphviz](https://graphviz.org/) to generate figures that show the relationship between species and reactions in a CRN. [`Graph(::SpeciesData, ::RxData)`](@ref) has been extended to allow for easy graph plotting directly from a Kinetica CRN:

```@example results_analysis
g = Graph(sd, rd)
savegraph(g, "../assets/tutorials/results_analysis/default_graph.svg", "svg"); nothing # hide
```

![](../assets/tutorials/results_analysis/default_graph.svg)

Graphs generated this way are compatible with Catalyst's graph methods and can be saved to file using `Catalyst.savegraph`. 

!!! note "Catalyst.jl Reexports"
    While Kinetica doesn't usually reexport any structs or methods from other packages, we make an exception for `Catalyst.Graph` and `Catalyst.savegraph` as we implement methods that can coexist and be interchanged with those in Catalyst. This allows [`Graph(::SpeciesData, ::RxData)`](@ref) to be called as above, but also allows `savegraph(::Graph, fname, fmt)` to be called without ever calling `using Catalyst` or needing it explicitly in your Julia project.

In addition, Kinetica allows for passing [Graphviz attributes](https://graphviz.org/doc/info/attrs.html) directly through to the plotter for extra customisability. These take the form of the `graph_attrs`, `species_attrs`, `rxn_attrs` and `edge_attrs` keyword arguments to [`Graph(::SpeciesData, ::RxData)`](@ref), which can be supplied with `Dict{Symbol, String}`s to modify parameters of the overall graph, the species nodes and the reaction nodes, and the edges respectively:

```@example results_analysis
g = Graph(sd, rd;
    graph_attrs=Dict(
        :layout => "sfdp",
        :overlap => "prism",
        :overlap_scaling => "-8"
    ), species_attrs=Dict(
        :shape => "hexagon",
        :color => "aquamarine1"
    ), rxn_attrs=Dict(
        :shape => "box",
        :color => "coral"
    ), edge_attrs=Dict(
        :color => "grey59"
    ))
savegraph(g, "../assets/tutorials/results_analysis/modified_graph.svg", "svg"); nothing # hide
```

![](../assets/tutorials/results_analysis/modified_graph.svg)

The exploration level that each species and reaction was found in (see [Iterative CRN Exploration](@ref)) is additionally added as a Graphviz attribute to each respective node. While this is (probably) not useful within the context of Graphviz's layout engines, more advanced graphical interfaces such as [Gephi](https://gephi.org/) can be used to create level-driven layouts:

![](../assets/crn11.png)

!!! note "On Graphviz Installation"
    Graphviz is one of Kinetica's Python dependencies, and as such is always installed at the same time as Kinetica. This is not the case for Catalyst.jl, which either requires a user-installed version of Graphviz, or it uses the [Graphviz_jll](https://github.com/JuliaBinaryWrappers/Graphviz_jll.jl) package. As discussed in the section on the [Graphviz](@ref) dependency, we avoid the JLL due to it missing some key features which Kinetica's CRNs typically need.

    Kinetica always adds its Python dependencies to the end of the current `PATH`, so if you have your own installation of Graphviz that you'd like to use, just make sure it's in your `PATH` before you start Julia!

### Species Analysis

Kinetica makes use of many functions within both [RDKit](https://github.com/rdkit/rdkit) and [Open Babel](https://github.com/openbabel/openbabel) to assist with basic property prediction and conversion between geometry and SMILES. A comprehensive list of functions implemented are available in the [Open Babel](@ref) and [RDKit](@ref) API pages.

Kinetica supports cached calculation of properties such as species molecular weights and hard-sphere radii through Open Babel through the [`get_species_stats!`](@ref) function. This places these properties in the `SpeciesData.cache`, which is useful for storing per-species values during a calculation.

Kinetica also stores species geometries within `SpeciesData.xyz`. These are stored as [ExtXYZ.jl](https://github.com/libAtoms/ExtXYZ.jl) `frame`s and can be manipulated as such.

## Analysing the Kinetic Simulation

Most of the information regarding the results of a kinetic simulation are stored in `ODESolveOutput.sol`, which is a DifferentialEquations `ODESolution` (in most cases, some extensions may apply depending on the simulation method but this can always be interacted as if it is an `ODESolution`). As such, the [DiffEq documentation](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/) applies for directly analysing simulation results.

However, when plotting simulation results, it can often be useful to have information about the experimental conditions that the simulation took place under, as well as some of the other parameters used in the simulation which are not included within the `ODESolution`. This includes the names of each of the species, which are only referred to be an integer ID within Kinetica's solution methods. Kinetica therefore includes [Plots.jl](https://github.com/JuliaPlots/Plots.jl) [plot recipes](https://docs.juliaplots.org/stable/recipes/) which act on the [`ODESolveOutput`](@ref) object, as this has access to all of these values alongside the `ODESolution`. The plots created by these recipes can be modified in the same way as any other Plots.jl figure.

### Concentration-Time Plots

This is the 'default' plot type for [`ODESolveOutput`](@ref), with the call signature:

```julia
plot(::ODESolveOutput; tunit="s", label_above=0.1, ignore_species=nothing, ignore_below=nothing)
```

where
* `tunit` is the unit of time to display on the x-axis;
* `label_above` is a concentration above which species should contribute to the plot's legend and be plotted in a colour (species is grey otherwise);
* `ignore_species` can be a `Vector` of SMILES strings representing species which should not be plotted (`nothing` if this is not desired);
* `ignore_below` can be a concentration threshold where species with a maximum concentration below are not plotted (`nothing` if this is not desired).

For example, plotting the same results as in [Getting Started](@ref), but excluding the radical species from the plot is simple:

```@setup results_analysis_2
using Kinetica

res = load_output("../my_CRN_out/direct_network_final.bson")
sd, rd = res.sd, res.rd
```

```@example results_analysis_2
using Plots

is_radical(smi) = ('[' in smi) && !(smi == "[H][H]")
radical_species = [spec for spec in [sd.toStr[i] for i in 1:sd.n] if is_radical(spec)]

plot(res; label_above=0.01, ignore_species=radical_species, linewidth=5)
savefig("../assets/tutorials/results_analysis/kinetics_plot.svg"); nothing # hide
```

![](../assets/tutorials/results_analysis/kinetics_plot.svg)

!!! note "Using Plots.jl Arguments"
    Notice how we've used the `linewidth` keyword argument here, but this isn't defined explicitly in this plot recipe. This is because Plots.jl allows for modifying figures, both within the original `plot` call and in subsequent `plot!` calls - we could've similarly done `plot!(linewidth=5)` after our original `plot` call to achieve the same goal.

### Conditions Plots

Kinetica defines a separate plot recipe for variable condition profiles, `conditionsplot`:

```julia
conditionsplot(::Union{::ODESolveOutput, ConditionSet}, ::Symbol; tunit="s")
```

This recipe can take either an [`ODESolveOutput`](@ref) or a [`ConditionSet`](@ref), but it has a second mandatory argument - a `Symbol` representing the condition profile to plot. For example, if a variable temperature profile is required (and it exists in the [`ConditionSet`](@ref) being plotted), then this is usually available under `:T`. Again, `tunit` is the unit of time to display on the x-axis.

```@example results_analysis_2
# Again, this is compatible with Plots.jl keyword arguments.
conditionsplot(res, :T; linecolor=:red)
savefig("../assets/tutorials/results_analysis/temperature_plot.svg"); nothing # hide
```

![](../assets/tutorials/results_analysis/temperature_plot.svg)

### Final Concentration Plots

It can often be useful to examine the final species concentrations at the end of a kinetic simulation, as these often correspond to experimental observables. Kinetica implements a specialised bar chart recipe for this, `finalconcplot`:

```julia
finalconcplot(::ODESolveOutput; 
              quantity=:conc, n_top=10, highlight_radicals=false, 
              ignore_species=nothing, xscale=:identity)
```

where
* `quantity` is a `Symbol`, either `:conc` or `:percent`. The former plots species concentrations, while the latter plots percentage concentrations with respect to the overall reaction mixture;
* `n_top` controls the number of species to include in the plot. These are sorted in decreasing order, with species outside `n_top` going into their own 'Others' bar at the end. For example, with `n_top=10`, the 10 species with the highest final concentrations will be plotted, followed by a bar for the combined final concentration of all other species;
* `highlight_radicals` colours radical species in red when `true`, as these are often not desired at the end of simulations;
* `ignore_species` accepts a `Vector` of SMILES strings for species to exclude from the plot. This can be useful when inert species are present with concentrations that are unimportant to the final state of the CRN (`nothing` if this is not desired);
* `xscale` mimics Plots.jl's `xscale` option for setting the scaling of the x-axis. In this case, only `:identity` (default) and `:log10` are supported, as additional calculations are required to obtain the correct x-axis limits for correct bar chart rendering.

```@example results_analysis_2
finalconcplot(res; quantity=:percent, n_top=5, highlight_radicals=true)
savefig("../assets/tutorials/results_analysis/concs_plot.svg"); nothing # hide
```

![](../assets/tutorials/results_analysis/concs_plot.svg)
