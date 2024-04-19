# [CRN Representation](@id crn_representation_page)

As noted in many of the main tutorials, Kinetica represents CRNs through two structs:

* [`SpeciesData`](@ref), which holds information on the species within a CRN including their SMILES representations, integer IDs and their geometries (stored as [ExtXYZ.jl](https://github.com/libAtoms/ExtXYZ.jl) `frame::Dict{String, Any}`).
* [`RxData`](@ref), which holds information about the reactions within a CRN including the integer IDs of the species within each reaction's reactants and products, the stoichiometries of these species, the atom-mapped reaction SMILES of each reaction, and unique hashes that can ensure reactions are not duplicated, irrespective of species ordering.

These are split into two rather than making a single unified CRN struct because:

* Situations often arise where only data regarding species or reactions are needed.
* This allows dedicated methods of base Julia functions such as `push!` and `splice!` to be easily added for adding/removing species and reactions.
* A unified CRN struct could easily become quite bloated.

## Representing Species (`SpeciesData`)

At the core of [`SpeciesData`](@ref) is a bidirectional mapping between species IDs and their respective SMILES representations within its `toStr` and `toInt` fields (which are just `Dict`s mapping ID to SMILES and SMILES to ID respectively). This mapping facilitates going back and forth between a human-understandable format and one that Kinetica can use when automatically assembling systems of ODEs for kinetic simulations.

!!! note "Species Within ODESystems"
    When a CRN gets converted to a ModelingToolkit `ODESystem`, every species needs to have a symbolic representation that can be mapped to its concentration during kinetic simulations. While it may be natural to assume we'd simply use the species' SMILES representations here, these contain special characters that cannot be present in `Symbol`s. To also facilitate simple scaling to any number of species, species are automatically mapped to time-dependent symbolic variables `spec(t)[i]`, where `i` becomes an integer ID within the range `1:SpeciesData.n`. This is one of the core reasons why this SMILES to ID mapping is required.

An instance of [`SpeciesData`](@ref) can be constructed in one of three ways:

* From an array of SMILES strings and a corresponding array of geometries (ExtXYZ `frame`s).
* From a single XYZ file containing one or more species.
* From a CDE mechanism exploration directory with [`import_mechanism`](@ref) (see later).
* As an empty struct (correctly initialised, but with no species), with [`init_network`](@ref).

Once created, more species can be added, again from SMILES and their respective geometries, or from other XYZ files, using various methods of [`push!(::SpeciesData, ...)`](@ref push!(::SpeciesData, ::String, ::Dict{String, Any})). These methods will all add species irrespective of whether they are already present. To only add new unique species, instead use the methods of [`push_unique!(::SpeciesData, ...)`](@ref push_unique!).

### Fields

Aside from the aforementioned `toStr` and `toInt` fields, [`SpeciesData`](@ref) structs can hold a wealth of information:

* `SpeciesData.n` holds the total number of species that have been added.
* `SpeciesData.xyz` holds a `Vector` of species geometries, indexed by species ID, in ExtXYZ.jl `frame` format. These store the cartesian coordinates of all atoms within a species, but can be used to store additional information related to a species.
* `SpeciesData.cache` holds an empty `Dict`, which can be used to store temporary information related to a species. For example, caling [`get_species_stats!`] populates the cache with `Dict`s that map species IDs to properties including molecular weight and hard-sphere radius.

### Use of SMILES

While [`SpeciesData`](@ref) admits any string-form identifier for species and uses this identifier as its method for testing if two species are identical, Kinetica as a whole currently uses SMILES as its only form of species identification. This is because SMILES is human-readable, and the vast majority of gas- and liquid-phase systems can be described using combinations of SMILES strings.

!!! note "Solid-Phase Systems"
    Support for periodic solid-phase systems in Kinetica is not currently planned and would require a different identifier that supports representation of periodic systems. However, with such an identifier in hand, Kinetica could be modified for this purpose.

Using SMILES as an identifier does come with a significant downside though - as a 2D descriptor of molecular structure, SMILES is not capable of distinguising between different conformers of molecular species. It is however capable of representing unique *cis*-*trans*/E-Z isomers and enantiomers. 

Since Kinetica is primarily made for investigating CRNs for long-timescale kinetic processes, we therefore view this as less of a downside of SMILES, and more of a necessary coarse-graining step. While molecular species usually require significant energy to move between the types of stereoisomer that SMILES *can* represent (e.g. temporary breaking of a ``\pi``-bond to move between *cis* and *trans* isomers around a C-C double bond), the purely conformational changes that SMILES *cannot* represent typically come with very low energy barriers. As such, over long timescales, contributions from any but the most stable conformation of a given species should be almost negligible, as any conformational changes should occur incredibly quickly.

Representing species as SMILES and assuming reactions only occur between the lowest energy conformer of each species therefore eliminates the need to consider these conformational changes. However, care should be taken when calculating reaction rate constants to always use low-energy conformations of species.

## Representing Reactions (`RxData`)

[`RxData`](@ref) maintains its space efficiency by only encoding reactions using species' integer IDs. They are therefore necessarily constructed in conjunction with an existing [`SpeciesData`](@ref), which **must** already contain the species that are in these reactions. They are therefore only constructed:

* As an empty struct, with [`init_network`](@ref).
* From arrays of reactant/product SMILES and reactant/product XYZ systems, alongside a [`SpeciesData`](@ref).

To avoid the hassle of importing reactions from files, adding their species to a [`SpeciesData`](@ref) and then adding their reactions to an [`RxData`](@ref), Kinetica alos provides shortcuts in the form of [`import_mechanism`](@ref) for starting new CRNs and [`import_mechanism!`](@ref) for extending existing ones. These will be covered in the [next section below](@ref "CRN File Structure").

To extend existing [`RxData`](@ref), [`push!(::RxData, ::SpeciesData, ...)`](@ref push!(::RxData, ::SpeciesData, reacs, prods, rsys, psys, dH)) should be used. This has an optional argument `unique_rxns=true` which can be used to ensure only unique reactions are added to a CRN.

To remove existing reactions from an [`RxData`](@ref), [`splice!(::RxData, ::Vector{Int})`](@ref) can be used. This allows users to specify a `Vector` of integer reaction IDs to be specified, and will remove these reactions from the CRN while shifting the remaining reactions back along each of the fields of [`RxData`](@ref) to form reactions with new IDs (much like how Julia's default `splice!` works on regular arrays).

### Fields

Instances of [`RxData`](@ref) contain the following fields:

* `RxData.nr` counts the number of reactions currently in the struct.
* `RxData.mapped_rxns` holds the atom-mapped reaction SMILES representations of each reaction in the struct. These can be useful for generating correctly atom-indexed reactant/product systems (see [Analysing the CRN](@ref) for an example).
* `id_reacs` and `id_prods` hold `Vector`s of the unique reactant/product IDs in each reaction. For example, if the CRN's [`SpeciesData`](@ref) dictated that species `"CC"` had ID `1` and `"[CH3]"` had ID `2`, and reaction `5` was `CC -> 2 [CH3]` (heterolytic cleavage of ethane), then `id_reacs[5] == ["CC"]` and `id_prods[5] == ["[CH3]"]`.
* `stoic_reacs` and `stoic_prods` hold `Vector`s of the stoichiometry of each species in the reactions in `id_reacs` and `id_prods` respectively. This is why in the above example, `id_prods[5]` does not contain two `"[CH3]"`s - because `stoic_reacs[5] == [1]` and `stoic_prods[5] == [2]`.
* `dH` holds each reaction's enthalpy of reaction, ``\Delta H_r``, as calculated by the level of electronic structure theory used within CDE during reaction generation.
* `rhash` holds a unique reaction hash for every reaction in the struct. This is calculated by taking the `sha256` of the concatenation of the sorted reactant and product SMILES, and is used for checking whether new reactions are duplicates of any others already in the struct.

### New Reaction Criteria

When adding new reactions to [`RxData`](@ref), Kinetica performs a few checks. Reactions are discarded if they:

* Only contain conformational changes. Species are internally represented by SMILES, so any changes not representable by SMILES (e.g. single bond rotations) lead to species being classed as the same. Note that E/Z isomers and enantiomers are distinct from one another in SMILES.
* Exceed the maximum molecularity. By default this is 2, so only unimolecular and bimolecular reactions are allowed to enter a CRN. This is controllable by the `max_molecularity` argument when constructing or `push!`ing to [`RxData`](@ref).
* Already exist in the [`RxData`](@ref) they are being added to, provided `unique_rxns=true`. This is checked by comparing a reaction's unique hash to those currently in the CRN.

## CRN File Structure

Chemical reactions are generated in Kinetica using [CDE](https://github.com/HabershonLab/cde). When calling CDE, XYZ files for each of the generated reactions are output into a given working directory. During CRN exploration, Kinetica automatically lays out a directory tree for these reactions to be explored in based on a main 'head' directory `rdir_head`, as well as the current `level` and `subspace` of exploration (see [Iterative CRN Exploration](@ref)) and the current iteration of exploration within that subspace `rcount`. At any given time during CRN exploration (expect in between levels, when kinetic simulations are being performed), Kinetica will be running CDE at

```julia
"$(rdir_head)/level_$(lpad(level, 3, "0"))/subspace_$(lpad(subspace, 3, "0"))/reac_$(lpad(rcount, 5, "0"))"
```

This location is held and accumulated within an instance of [`Kinetica.ExploreLoc`](@ref) during CRN explorations, which implements a method [`pathof(::Kinetica.ExploreLoc)`](@ref) to quickly create this path.

Once a given CDE iteration is complete and there are new reactions at the current `rcount`, Kinetica reads these reactions in and adds them to the current CRN. It does this by calling [`import_mechanism!`](@ref), which does the following:

* Reads in all final reaction XYZs from the current exploration directory and processes them into separable reactant/product SMILES and ExtXYZ `frame`s using [`ingest_cde_run`](@ref).
* Adds all newly discovered species to the current [`SpeciesData`](@ref) using [`push_unique!(::SpeciesData)`](@ref push_unique!).
* Adds all newly discovered reactions to the current [`RxData`](@ref) using [`push!(::RxData, ::SpeciesData, ...; unique_rxns=true)`](@ref push!(::RxData, ::SpeciesData, reacs, prods, rsys, psys, dH)).

!!! note "Reverse Reactions"
    When [`ingest_cde_run`](@ref) is called to import all reactions in a CDE run into Kinetica, it automatically adds their reverse reactions too! In this way, we preserve the principle of detailed balance within our CRN.

### Checkpoints and Status Files

Once exploration within a subspace is complete, Kinetica writes a small file called `isconv` to the subspace's directory. This acts as a checkpoint signal that tells Kinetica that this subspace is converged (completed) if it ever needs to re-import this CRN from its file tree. This can occur when CRN exploration is ended early, for example due to an error or a lack of walltime, and the Julia process is ended.

If the same exploration script is run again and the directory structure is unchanged, Kinetica will be able to pick up on the latest `isconv` file and know that this was the last completed subspace with [`Kinetica.find_current_loc`](@ref). It can then import the currently explored CRN with [`import_network`](@ref) and continue from the next subspace.

Similarly, during iterative CRN explorations, Kinetica leaves behind files containing information about the seed species that should form the subspaces in each level. These `seeds.in` files are placed in each level's directory upon level creation, and similarly act as important markers for determination of exploration location within [`Kinetica.find_current_loc`](@ref). However, they are mainly used so that Kinetica can know which seed species have had their self-reactions fully explored in previous levels ([`Kinetica.load_past_seeds`](@ref)), and which seed species are responsible for the subspaces in the current level ([`Kinetica.load_current_seeds`](@ref)).

Once a level's exploration is complete and Kinetica has finished a kinetic simulation to determine the next seed species, it calls [`Kinetica.identify_next_seeds`](@ref). This uses concentration criteria defined in [`IterativeExplore`](@ref) to choose species with kinetic relevance to the overall CRN. Upon identification of these seeds, Kinetica writes a `seeds.out` file to the current exploration's `savedir` (provided `savedir` is given as a keyword argument when calling [`explore_network`](@ref)). This file is not currently used for restarting explorations, but can be used to monitor the progress of the current CRN exploration.