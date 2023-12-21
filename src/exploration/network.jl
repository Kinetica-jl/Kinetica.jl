"""
Bidirectional String-Int dictionary for chemical species.

Contains fields for:
* SMILES string -> integer ID dictionary (`toInt`)
* Integer ID -> SMILES string dictionary (`toStr`)
* Number of species (`n`)
* ExtXYZ structures of species (`xyz`)
* Dictionary of per-species cached values (`cache`)
"""
mutable struct SpeciesData{iType}
    toInt::Dict{String, iType}
    toStr::Dict{iType, String}
    n::iType
    xyz::Dict{iType, Dict{String, Any}}
    cache::Dict{Any, Any}
end

"""
    sd = SpeciesData(smi_list, xyz_list)

Outer constructor method for `SpeciesData`, allowing for construction
from a list of SMILES strings and ExtXYZ structures.
"""
function SpeciesData(smi_list, xyz_list; unique_species=true)
    n = length(smi_list)
    if n == 0
        return SpeciesData(Dict(), Dict(), 0, Dict(), Dict())
    end

    if unique_species
        unique_smi_list = []
        unique_xyz_list = []
        for i in 1:n
            if !(smi_list[i] in unique_smi_list)
                push!(unique_smi_list, smi_list[i])
                push!(unique_xyz_list, xyz_list[i])
            end
        end
        return SpeciesData(
            Dict(smi => i for (i, smi) in enumerate(unique_smi_list)),
            Dict(i => smi for (i, smi) in enumerate(unique_smi_list)),
            length(unique_smi_list),
            Dict(i => x for (i, x) in enumerate(unique_xyz_list)),
            Dict()
        )
    else
        return SpeciesData(
            Dict(smi => i for (i, smi) in enumerate(smi_list)),
            Dict(i => smi for (i, smi) in enumerate(smi_list)),
            n,
            Dict(i => x for (i, x) in enumerate(xyz_list)),
            Dict()
        )
    end
end

"""
    sd = SpeciesData(xyz_file[, fix_radicals])

Outer constructor method for `SpeciesData`, allowing for construction
from an XYZ file.
"""
function SpeciesData(xyz_file::String; fix_radicals=true)
    smi_list, xyz_list = ingest_xyz_system(xyz_file; fix_radicals)
    SpeciesData(smi_list, xyz_list)
end

"""
    push!(sd, smi, xyz)

Add a species to `SpeciesData`.

Does not account for `smi` already existing within `sd`. To
ensure no overlap, use `push_unique!`.
"""
function Base.push!(sd::SpeciesData, smi::String, xyz::Dict{String, Any})
    sd.n += 1
    sd.toInt[smi] = sd.n
    sd.toStr[sd.n] = smi
    sd.xyz[sd.n] = xyz
    return
end

"""
    push!(sd, xyz_file[, fix_radicals])

Add all species in `xyz_file` to `sd`.

Does not account for `smi` already existing within `sd`. To
ensure no overlap, use `push_unique!`.
"""
function Base.push!(sd::SpeciesData, xyz_file::String; fix_radicals=true)
    smi_list, xyz_list = ingest_xyz_system(xyz_file; fix_radicals)
    for (smi, xyz) in zip(smi_list, xyz_list)
        push!(sd, smi, xyz)
    end
    return
end

"""
    push!(sd, smis, xyzs)

Add an array of species to `SpeciesData`.

Does not account for `smi` already existing within `sd`. To
ensure no overlap, use `push_unique!`.
"""
function Base.push!(sd::SpeciesData, smis::Vector{String}, xyzs::Vector{Any})
    for (smi, xyz) in zip(smis, xyzs)
        push!(sd, smi, xyz)
    end
    return
end

"""
    push_unique!(sd, smi, xyz)

Add a species SMILES to a `SpeciesData`, as long as it does not already exist there.
"""
function push_unique!(sd::SpeciesData, smi::String, xyz::Dict{String, Any})
    if !(smi in keys(sd.toInt))
        push!(sd, smi, xyz)
    end
    return
end

"""
    push_unique!(sd, xyz_file[, fix_radicals])

Add species in `xyz_file` to `sd`, as long as they do not already exist there.
"""
function push_unique!(sd::SpeciesData, xyz_file::String; fix_radicals=true)
    smi_list, xyz_list = ingest_xyz_system(xyz_file; fix_radicals)
    for (smi, xyz) in zip(smi_list, xyz_list)
        if !(smi in keys(sd.toInt))
            push!(sd, smi, xyz)
        end
    end
    return
end

"""
    push_unique!(sd, smis, xyzs)

Add an array of species to `SpeciesData`, as long as each does not already exist there.
"""
function push_unique!(sd::SpeciesData, smis::Vector{String}, xyzs::Vector{Dict{String, Any}})
    for (smi, xyz) in zip(smis, xyzs)
        push_unique!(sd, smi, xyz)
    end
    return
end


"""
Data container for reactions.

Should be constructed alongside a `SpeciesData` for mapping
species IDs to SMILES strings.

Contains fields for:
* Number of reactions encountered (`nr`)
* Atom-mapped reaction SMILES for unambiguous linking of atom indices in reactants and products (`mapped_rxns`)
* Unique IDs of reactants for each reaction (`id_reacs`)
* Unique IDs of products for each reaction (`id_prods`)
* Stoichiometries of reactants for each reaction (`stoic_reacs`)
* Stoichiometries of products for each reaction (`stoic_prods`)
* Reaction enthalpies (`dH`)
* Reaction hashes, used for unique identification (`rhash`)
"""
mutable struct RxData{iType, fType}
    nr::iType
    mapped_rxns::Vector{String}
    id_reacs::Vector{Vector{iType}}
    id_prods::Vector{Vector{iType}}
    stoic_reacs::Vector{Vector{iType}}
    stoic_prods::Vector{Vector{iType}}
    dH::Vector{fType}
    rhash::Vector{Vector{UInt8}}
end

"""
    rd = RxData(sd, reacs, prods, rsys, psys, dH[, unique_rxns, max_molecularity])

Outer constructor method for `RxData`.

Creates a new `RxData` reaction data store from lists of
reactant and product SMILES by cross-referencing species
IDs with `sd`. `sd` should therefore already have all
species input here loaded in.

`reacs` and `prods` should be the raw SMILES
arrays from `ingest_cde_run()`, i.e. without any duplicate
species removed due to stoichiometry. This constructor will
determine stoichiometry and output the unique form of each
set of reactants/products.

By default, adds a maximum of 1 of each reaction type when
`unique_rxns = true`, and only admits reactions with a maximum
molecularity of 2 (i.e. bimolecular reactions).
"""
function RxData(sd::SpeciesData{iType}, 
        reacs::Vector{Vector{String}}, prods::Vector{Vector{String}},
        rsys::Vector{Dict{String, Any}}, psys::Vector{Dict{String, Any}},
        dH::Vector{fType}; 
        unique_rxns=true, max_molecularity=2) where {iType, fType <: AbstractFloat}

    rxns_final = []
    id_reacs_final = []
    id_prods_final = []
    stoic_reacs_final = []
    stoic_prods_final = []
    dH_final = []
    hashes_final = []
    nr = 0

    invcounter = 0
    dcounter = 0
    for i in 1:length(reacs)
        # Obtain accumulators so that we have unique sets of species with counts.
        reac_counter = counter(reacs[i])
        prod_counter = counter(prods[i])

        # Check for purely conformational changes, which are invalid.
        if issetequal(reac_counter, prod_counter)
            invcounter += 1
            continue
        end

        # Check there are no reactions exceeding max_molecularity (forward or backward).
        if length(reac_counter) > max_molecularity || length(prod_counter) > max_molecularity || 
            sum(values(reac_counter)) > max_molecularity || sum(values(prod_counter)) > max_molecularity
            invcounter += 1
            continue
        end

        # Obtain reaction hash.
        all_reacs = sort(reduce(vcat, [[key for _ in 1:reac_counter[key]] for key in keys(reac_counter)]))
        all_prods = sort(reduce(vcat, [[key for _ in 1:prod_counter[key]] for key in keys(prod_counter)]))
        rhash = stable_hash(vcat(all_reacs, all_prods))

        # Add reaction to arrays if it is unique.
        if !unique_rxns || !(rhash in hashes_final)
            # Construct atom-mapped reaction SMILES from original geometries.
            mapped_reacs = atom_map_smiles(rsys[i], join(all_reacs, "."))
            mapped_prods = atom_map_smiles(psys[i], join(all_prods, "."))
            mapped_rxn = join([mapped_reacs, mapped_prods], ">>")

            # Get unique lists so stoichiometry can multiply each species in a reaction.
            unique_reacs = unique(all_reacs)
            unique_prods = unique(all_prods)

            nr += 1
            push!(rxns_final, mapped_rxn)
            push!(id_reacs_final, [sd.toInt[reac] for reac in unique_reacs])
            push!(id_prods_final, [sd.toInt[prod] for prod in unique_prods])
            push!(stoic_reacs_final, [reac_counter[reac] for reac in unique_reacs])
            push!(stoic_prods_final, [prod_counter[prod] for prod in unique_prods])
            push!(dH_final, dH[i])
            push!(hashes_final, rhash)
        else
            dcounter += 1
        end
    end

    @debug " - $dcounter duplicate and $invcounter invalid reactions found."

    return RxData{iType, fType}(nr, rxns_final, id_reacs_final, id_prods_final,
            stoic_reacs_final, stoic_prods_final, dH_final, hashes_final)
end

"""
    push!(rd, sd, reacs, prods, rsys, psys, dH[, unique_rxns, max_molecularity])
"""
function Base.push!(rd::RxData{iType, fType}, sd::SpeciesData,
        reacs::Vector{Vector{String}}, prods::Vector{Vector{String}}, 
        rsys::Vector{Dict{String, Any}}, psys::Vector{Dict{String, Any}},
        dH::Vector{fType}; 
        unique_rxns=true, max_molecularity=2) where {iType, fType <: AbstractFloat}

    invcounter = 0
    dcounter = 0
    for i in 1:length(reacs)
        # Obtain accumulators so that we have unique sets of species with counts.
        reac_counter = counter(reacs[i])
        prod_counter = counter(prods[i])

        # Check for purely conformational changes, which are invalid.
        if issetequal(reac_counter, prod_counter)
            invcounter += 1
            continue
        end

        # Check there are no reactions exceeding max_molecularity (forward or backward).
        if length(reac_counter) > max_molecularity || length(prod_counter) > max_molecularity || 
            sum(values(reac_counter)) > max_molecularity || sum(values(prod_counter)) > max_molecularity
            invcounter += 1
            continue
        end

        # Obtain reaction hash.
        all_reacs = sort(reduce(vcat, [[key for _ in 1:reac_counter[key]] for key in keys(reac_counter)]))
        all_prods = sort(reduce(vcat, [[key for _ in 1:prod_counter[key]] for key in keys(prod_counter)]))
        rhash = stable_hash(vcat(all_reacs, all_prods))

        # Add reaction to arrays if it is unique.
        if !unique_rxns || !(rhash in rd.rhash)
            # Construct atom-mapped reaction SMILES from original geometries.
            mapped_reacs = atom_map_smiles(rsys[i], join(all_reacs, "."))
            mapped_prods = atom_map_smiles(psys[i], join(all_prods, "."))
            mapped_rxn = join([mapped_reacs, mapped_prods], ">>")

            # Get unique lists so stoichiometry can multiply each species in a reaction.
            unique_reacs = unique(all_reacs)
            unique_prods = unique(all_prods)

            rd.nr += 1
            push!(rd.mapped_rxns, mapped_rxn)
            push!(rd.id_reacs, [sd.toInt[reac] for reac in unique_reacs])
            push!(rd.id_prods, [sd.toInt[prod] for prod in unique_prods])
            push!(rd.stoic_reacs, [reac_counter[reac] for reac in unique_reacs])
            push!(rd.stoic_prods, [prod_counter[prod] for prod in unique_prods])
            push!(rd.dH, dH[i])
            push!(rd.rhash, rhash)
            
        else
            dcounter += 1
        end
    end

    @debug " - $dcounter duplicate and $invcounter invalid reactions found."

    return
end


"""
    sd, rd = init_network([iType, fType])

Initialises an empty reaction network.
"""
function init_network(; iType=Int64, fType=Float64)
    sd = SpeciesData{iType}(
        Dict{String, iType}(), Dict{iType, String}(),
        0, 
        Dict{iType, Dict{String, Any}}(), Dict()
    )
    rd = RxData{iType, fType}(
        0, String[],
        Vector{iType}[], Vector{iType}[], 
        Vector{iType}[], Vector{iType}[], 
        fType[], Vector{UInt8}[]
    )

    return sd, rd
end


"""
    splice!(rd, rids)

Removes reactions at indeces `rids` from `rd`.
"""
function Base.splice!(rd::RxData, rids::Vector{Int})
    if current_logger().min_level <= Debug
        for rid in rids
            @debug "Removing reaction $rid from network:"
        end
    end

    if length(rids) > 0
        for f in fieldnames(RxData)
            if f != :nr
                splice!(getfield(rd, f), rids)
            end
        end
        rd.nr -= length(rids)
    end
end


"""
    rxn_str = format_rxn(sd, rd, rid)

Nicely formats a string describing the reaction at `rid`.
"""
function format_rxn(sd::SpeciesData, rd::RxData, rid::Int)
    reacs = [sd.toStr[sid] for sid in rd.id_reacs[rid]]
    prods = [sd.toStr[sid] for sid in rd.id_prods[rid]]
    reac_strs = [n > 1 ? "$n $spec" : spec for (n, spec) in zip(rd.stoic_reacs[rid], reacs)]
    prod_strs = [n > 1 ? "$n $spec" : spec for (n, spec) in zip(rd.stoic_prods[rid], prods)]
    rxn_str = join([join(reac_strs, " + "), join(prod_strs, " + ")], " --> ")
    return rxn_str
end