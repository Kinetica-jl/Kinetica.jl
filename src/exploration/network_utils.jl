"""
    get_rhash(sd::SpeciesData, rd::RxData, rid)

Returns the reaction hash for the reaction at `rid` in `rd`.
"""
function get_rhash(sd::SpeciesData, rd::RxData, rid)
    reacs = String[]
    for (i, sid) in enumerate(rd.id_reacs[rid])
        for _ in 1:rd.stoic_reacs[rid][i] 
            push!(reacs, sd.toStr[sid])
        end
    end
    sort!(reacs)
    prods = String[]
    for (i, sid) in enumerate(rd.id_prods[rid])
        for _ in 1:rd.stoic_prods[rid][i] 
            push!(prods, sd.toStr[sid])
        end
    end
    sort!(prods)

    return stable_hash(vcat(reacs, prods); version=4)
end


"""
    get_reverse_rhash(sd::SpeciesData, rd::RxData, rid)

Returns the reverse reaction hash for the reaction at `rid` in `rd`.

Useful when needing to identify if a reverse reaction
is already in a CRN without having to look through many
species permutations.
"""
function get_reverse_rhash(sd::SpeciesData, rd::RxData, rid)
    reacs = String[]
    for (i, sid) in enumerate(rd.id_reacs[rid])
        for _ in 1:rd.stoic_reacs[rid][i] 
            push!(reacs, sd.toStr[sid])
        end
    end
    sort!(reacs)
    prods = String[]
    for (i, sid) in enumerate(rd.id_prods[rid])
        for _ in 1:rd.stoic_prods[rid][i] 
            push!(prods, sd.toStr[sid])
        end
    end
    sort!(prods)

    forw_rhash = stable_hash(vcat(reacs, prods); version=4)
    @assert rd.rhash[rid] == forw_rhash # Checks StableHashTraits hasn't broken anything.
    rev_rhash = stable_hash(vcat(prods, reacs); version=4)
    return rev_rhash
end


"""
    get_species_stats!(sd::SpeciesData[, refresh=false])

Gets statistics about the species in `sd`.

Calculates values for the following useful statistics, which
are placed in the species data cache at `sd.cache` with a key
for each property:

* Average COM-atom radius of species (`:radii`)
* Molecular weight of species (`:weights`)

Only calculates statistics for species without defined values,
unless `refresh=true` which causes all statistics for every
species to be updated.
"""
function get_species_stats!(sd::SpeciesData{iType}; refresh::Bool=false) where {iType}
    # Create cache dicts if they don't already exist.
    properties = [:radii, :weights]
    property_types = [Float64, Float64]
    for (prop, pType) in zip(properties, property_types)
        if !(prop in keys(sd.cache))
            sd.cache[prop] = Dict{iType, pType}()
        end
    end

    for i in 1:sd.n
        if refresh || !all([i in keys(sd.cache[prop]) for prop in properties])
            atoms = frame_to_atoms(sd.xyz[i])
            sd.cache[:weights][i] = pyconvert(Float64, sum(atoms.get_masses()))
            na = sd.xyz[i]["N_atoms"]
            if na == 1
                sd.cache[:radii][i] = pyconvert(Float64, ase.data.vdw_radii[atoms.get_atomic_numbers()[1]])
            else
                sd.cache[:radii][i] = calc_average_molecular_radius(atoms)
            end
        end
    end
end


"""
    calc_average_molecular_radius(atoms)

Calculates average COM-atom radius of a given ASE Atoms object.

Calculates the radial distance from the center of mass (COM)
to each atom, then takes the mean of these values to obtain
an average molecular radius. Adds on a correction based on the
average van der Walls radius of the atoms to emulate the outer
edge of the molecule.
"""
function calc_average_molecular_radius(atoms)
    coords = pyconvert(Matrix{Float64}, atoms.get_positions() - atoms.get_center_of_mass())
    r_avg = mean(norm.(eachrow(coords)))
    avg_vdw = mean(pyconvert(Vector{Float64}, [ase.data.vdw_radii[n] for n in atoms.get_atomic_numbers()]))
    radius = r_avg + avg_vdw
    return radius
end


"""
    format_rxn(sd::SpeciesData, rd::RxData, rid::Int[, display_level=false])

Nicely formats a string describing the reaction at `rid`.

If the keyword argument `display_level` is `true`, additionally
annotates the reaction with the level in which it was discovered.
"""
function format_rxn(sd::SpeciesData, rd::RxData, rid::Int; display_level=false)
    reacs = [sd.toStr[sid] for sid in rd.id_reacs[rid]]
    prods = [sd.toStr[sid] for sid in rd.id_prods[rid]]
    reac_strs = [n > 1 ? "$n $spec" : spec for (n, spec) in zip(rd.stoic_reacs[rid], reacs)]
    prod_strs = [n > 1 ? "$n $spec" : spec for (n, spec) in zip(rd.stoic_prods[rid], prods)]
    rxn_str = join([join(reac_strs, " + "), join(prod_strs, " + ")], " --> ")
    if display_level
        rxn_str = "L$(rd.level_found[rid]): " * rxn_str
    end
    return rxn_str
end

"""
    print_rxn(sd::SpeciesData, rd::RxData, rid::Int[, display_level=false])

Prints the reaction at ID `rid` with SMILES names for species.

If the keyword argument `display_level` is `true`, additionally
annotates the reaction with the level in which it was discovered.
"""
function print_rxn(sd::SpeciesData, rd::RxData, rid::Int; display_level=false)
    println(format_rxn(sd, rd, rid; display_level=display_level))
end