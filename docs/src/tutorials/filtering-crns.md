# Filtering CRNs

Sometimes, certain chemical reactions are known to be unwanted in a CRN before any kinetic simulations are run. For example, it may be desirable to place a limit on the size of chemical species within a CRN by removing any reactions which create species above this size.

Kinetica enables this reaction filtering through an easily extensible system that can scale to as many arbitrary filters as desired. This is done by constructing a [`RxFilter`](@ref). This is essentially a `Vector` of functions, each of which takes the current CRN in the form `(SpeciesData, RxData)` and returns a boolean mask over the reactions in the given [`RxData`](@ref).

To demonstrate, let's construct an [`RxFilter`](@ref) which removes all reactions involving double bonds from a CRN. We'll use the CRN we created in [Getting Started](@ref) by loading the results file and extracting the [`SpeciesData`](@ref) and [`RxData`](@ref):

```@example filters
using Kinetica
res = load_output("../my_CRN_out/direct_network_final.bson")
sd, rd = res.sd, res.rd
nothing # hide
```

If we inspect the reactions in this CRN, we can see that quite a few of them involve double bonds (represented in SMILES by an equals sign `=` between two elements):

```@example filters
for i in 1:rd.nr
    print_rxn(sd, rd, i)
end
nothing # hide
```

To create our [`RxFilter`](@ref), we must first write a function that identifies reactions which involve double bonds, then writes a boolean mask which indicates that these reactions should be removed. Remember that this function must take the current CRN in the form `(SpeciesData, RxData)`:

```@example filters
function db_filter(sd, rd)
    # Create a mask where no reactions are removed.
    mask = [false for _ in 1:rd.nr]

    # Iterate through products of reactions.
    for (i, rx_reac_ids) in enumerate(rd.id_reacs)
        for reac_species_id in rx_reac_ids
            # Set mask at this reaction to true if there are double bonds.
            # We can determine this easily by the presence of '=' in a species' SMILES.
            if '=' in sd.toStr[reac_species_id]
                mask[i] = true
                break
            end
        end
    end
    
    # Do the same for products of reactions.
    for (i, rx_prod_ids) in enumerate(rd.id_prods)
        # Skip already masked reactions.
        if mask[i] continue end
        for prod_species_id in rx_prod_ids
            if '=' in sd.toStr[prod_species_id]
                mask[i] = true
                break
            end
        end
    end

    return mask
end
```

We can turn this into a [`RxFilter`](@ref) like so:

```@example filters
filter = RxFilter([db_filter])
```

!!! note "Inverting Filters"
    We have created a filter that will remove reactions where the mask is `true`, however sometimes it may be more convenient to write filters that only *keep* reactions where the mask is `true`. This is supported by passing the keyword argument `keep_filtered=true` (`false` by default) to [`RxFilter`](@ref) during its construction.

Filters are applied automatically when running a kinetic simulation, as long as they are provided when constructing the intermediate [`StaticODESolve`](@ref)/[`VariableODESolve`](@ref) object:

```julia
# conditions = ConditionSet(...)
# pars = ODESimulationParams(...)
# calculator = ...
solvemethod = VariableODESolve(pars, conditions, calculator, filter)
# res = solve_network(solvemethod, sd, rd)
```

Under the hood, this is just calling a simple procedure: `filter` is queried to construct a combined mask from all filter functions given the current CRN, and [`splice!(::RxData, ::Vector{Int})`](@ref) is called to remove the masked reactions from the CRN. We can replicate this to see how the CRN is modified by our filter:

```@example filters
mask = get_filter_mask(filter, sd, rd)
# Create a copy of RxData to avoid modifying it.
filtered_rd = deepcopy(rd)
splice!(filtered_rd, findall(mask))

for i in 1:filtered_rd.nr
    print_rxn(sd, filtered_rd, i)
end
nothing # hide
```

!!! note "Using findall(mask)"
    Note that we can't just pass the `mask` directly to the `splice!` function - instead, we must pass it through `findall`. This converts the `BitVector` `mask` into a `Vector{Int}` containing the indices of reactions where the mask is `true`, which is compatible with `splice!`.

We've now successfully masked out all the reactions that involve double bonds from our CRN! However, let's say we also want to remove all reactions that consume methane (`C` in SMILES). We could extend our `db_filter` function above, but it would start to become quite large. Instead, we can create a new function for this new filter and let the [`RxFilter`](@ref) do the work of combining the resulting masks together:

```@example filters
function methane_filter(sd, rd)
    mask = [false for _ in 1:rd.nr]
    for (i, rx_reac_ids) in enumerate(rd.id_reacs)
        for reac_species_id in rx_reac_ids
            if sd.toStr[reac_species_id] == "C"
                mask[i] = true
                break
            end
        end
    end
    return mask
end

# Create a new filter that uses both filter functions.
filter = RxFilter([db_filter, methane_filter])

# Create a new combined mask to act on the original CRN.
mask = get_filter_mask(filter, sd, rd)

# Filter the original CRN with both filters.
filtered_rd = deepcopy(rd)
splice!(filtered_rd, findall(mask))

for i in 1:filtered_rd.nr
    print_rxn(sd, filtered_rd, i)
end
nothing # hide
```

Notice that now we've also filtered out all reactions where methane (`C`) was a reactant! We can combine as many filter functions as we want in this way to remove arbitrary sets of reactions from a CRN as needed.