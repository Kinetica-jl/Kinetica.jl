mutable struct SpeciesData{iType}
    toInt::Dict{String, iType}
    toStr::Dict{iType, String}
    n::iType
    xyz::Dict{iType, Dict{String, Any}}
    level_found::Dict{iType, Int}
    cache::Dict{Any, Any}
end

"""
    SpeciesData(smi_list, xyz_list[, unique_species=true])
    SpeciesData(smi_list, xyz_list, level[, unique_species=true])
    SpeciesData(xyz_file[, unique_species=true, fix_radicals=true])
    SpeciesData(xyz_file, level[, unique_species=true, fix_radicals=true])

Bidirectional String-Int dictionary for chemical species.

Can either be constructed from an array of SMILES strings
and their corresponding ExtXYZ frames, or from a single
XYZ file with one or multiple species present. If
`unique_species=true`, will not include any duplicate
species if present. If `fix_radicals=true` in the XYZ
file loading case, will attempt to tidy up radical
SMILES with OBCR.

Both constuctors can optionally be given a `level` argument,
indicating the exploration level which a species (or set of
species) was first discovered. If this is not provided,
assumes species are entering the struct at level 1.

Contains fields for:
* SMILES string -> integer ID dictionary (`toInt`)
* Integer ID -> SMILES string dictionary (`toStr`)
* Number of species (`n`)
* ExtXYZ structures of species (`xyz`)
* Integer ID -> initial discovered level dictionary (`level_found`)
* Dictionary of per-species cached values (`cache`)
"""
function SpeciesData(smi_list, xyz_list, level; unique_species=true)
    n = length(smi_list)
    if n == 0
        return SpeciesData(Dict(), Dict(), 0, Dict(), Dict(), Dict())
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
            Dict(i => level for i in 1:length(unique_smi_list)),
            Dict()
        )
    else
        return SpeciesData(
            Dict(smi => i for (i, smi) in enumerate(smi_list)),
            Dict(i => smi for (i, smi) in enumerate(smi_list)),
            n,
            Dict(i => x for (i, x) in enumerate(xyz_list)),
            Dict(i => level for i in 1:n),
            Dict()
        )
    end
end
SpeciesData(smi_list, xyz_list; unique_species=true) = SpeciesData(smi_list, xyz_list, 1; unique_species=unique_species)

function SpeciesData(xyz_file::String, level::Int; unique_species=true, fix_radicals=true)
    smi_list, xyz_list = ingest_xyz_system(xyz_file; fix_radicals)
    SpeciesData(smi_list, xyz_list, level; unique_species=unique_species)
end
SpeciesData(xyz_file::String; unique_species=true, fix_radicals=true) = SpeciesData(xyz_file, 1; unique_species=unique_species, fix_radicals=fix_radicals)

"""
    push!(sd::SpeciesData, smi::String, xyz::Dict{String, Any})
    push!(sd::SpeciesData, smi::String, xyz::Dict{String, Any}, level::Int)

Add a species to `SpeciesData`.

Optionally takes a specified exploration level (defaults to 1
if not provided). Does not account for `smi` already existing
within `sd`. To ensure no overlap, use `push_unique!`.
"""
function Base.push!(sd::SpeciesData, smi::String, xyz::Dict{String, Any}, level::Int)
    sd.n += 1
    sd.toInt[smi] = sd.n
    sd.toStr[sd.n] = smi
    sd.xyz[sd.n] = xyz
    sd.level_found[sd.n] = level
    return
end
Base.push!(sd::SpeciesData, smi::String, xyz::Dict{String, Any}) = push!(sd, smi, xyz, 1)

"""
    push!(sd::SpeciesData, xyz_file::String[, fix_radicals=true])
    push!(sd::SpeciesData, xyz_file::String, level::Int[, fix_radicals=true])

Add all species in `xyz_file` to `sd`.

Optionally takes a specified exploration level (defaults to 1
if not provided). Does not account for `smi` already existing
within `sd`. To ensure no overlap, use `push_unique!`.
"""
function Base.push!(sd::SpeciesData, xyz_file::String, level::Int; fix_radicals=true)
    smi_list, xyz_list = ingest_xyz_system(xyz_file; fix_radicals)
    for (smi, xyz) in zip(smi_list, xyz_list)
        push!(sd, smi, xyz, level)
    end
    return
end
Base.push!(sd::SpeciesData, xyz_file::String; fix_radicals=true) = push!(sd, xyz_file, 1; fix_radicals=fix_radicals)

"""
    push!(sd::SpeciesData, smis::Vector{String}, xyzs::Vector{Any})
    push!(sd::SpeciesData, smis::Vector{String}, xyzs::Vector{Any}, level::Int)

Add an array of species to `SpeciesData`.

Optionally takes a specified exploration level (defaults to 1
if not provided). Does not account for `smi` already existing
within `sd`. To ensure no overlap, use `push_unique!`.
"""
function Base.push!(sd::SpeciesData, smis::Vector{String}, xyzs::Vector{Any}, level::Int)
    for (smi, xyz) in zip(smis, xyzs)
        push!(sd, smi, xyz, level)
    end
    return
end
Base.push!(sd::SpeciesData, smis::Vector{String}, xyzs::Vector{Any}) = push!(sd, smis, xyzs, 1)

"""
    push_unique!(sd::SpeciesData, smi::String, xyz::Dict{String, Any})
    push_unique!(sd::SpeciesData, smi::String, xyz::Dict{String, Any}, level::Int)

Add a species SMILES to a `SpeciesData`, as long as it does not already exist there.

Optionally takes a specified exploration level (defaults to 1
if not provided).
"""
function push_unique!(sd::SpeciesData, smi::String, xyz::Dict{String, Any}, level::Int)
    if !(smi in keys(sd.toInt))
        push!(sd, smi, xyz, level)
    end
    return
end
push_unique!(sd::SpeciesData, smi::String, xyz::Dict{String, Any}) = push_unique!(sd, smi, xyz, 1)

"""
    push_unique!(sd::SpeciesData, xyz_file::String[, fix_radicals=true])
    push_unique!(sd::SpeciesData, xyz_file::String, level::Int[, fix_radicals=true])

Add species in `xyz_file` to `sd`, as long as they do not already exist there.

Optionally takes a specified exploration level (defaults to 1
if not provided).
"""
function push_unique!(sd::SpeciesData, xyz_file::String, level::Int; fix_radicals=true)
    smi_list, xyz_list = ingest_xyz_system(xyz_file; fix_radicals)
    for (smi, xyz) in zip(smi_list, xyz_list)
        if !(smi in keys(sd.toInt))
            push!(sd, smi, xyz, level)
        end
    end
    return
end
push_unique!(sd::SpeciesData, xyz_file::String; fix_radicals=true) = push_unique!(sd, xyz_file, 1; fix_radicals=fix_radicals)

"""
    push_unique!(sd::SpeciesData, smis::String, xyzs::Vector{Dict{String, Any}})
    push_unique!(sd::SpeciesData, smis::String, xyzs::Vector{Dict{String, Any}}, level::Int)

Add an array of species to `SpeciesData`, as long as each does not already exist there.

Optionally takes a specified exploration level (defaults to 1
if not provided).
"""
function push_unique!(sd::SpeciesData, smis::Vector{String}, xyzs::Vector{Dict{String, Any}}, level::Int)
    for (smi, xyz) in zip(smis, xyzs)
        push_unique!(sd, smi, xyz, level)
    end
    return
end
push_unique!(sd::SpeciesData, smis::Vector{String}, xyzs::Vector{Dict{String, Any}}) = push_unique!(sd, smis, xyzs, 1)


mutable struct RxData{iType, fType}
    nr::iType
    mapped_rxns::Vector{String}
    id_reacs::Vector{Vector{iType}}
    id_prods::Vector{Vector{iType}}
    stoic_reacs::Vector{Vector{iType}}
    stoic_prods::Vector{Vector{iType}}
    dH::Vector{fType}
    rhash::Vector{Vector{UInt8}}
    level_found::Vector{Int}
end

"""
    rd = RxData(sd::SpeciesData, reacs::Vector{Vector{String}}, prods::Vector{Vector{String}}, 
                rsys::Vector{Dict{String, Any}}, psys::Vector{Dict{String, Any}}, 
                dH::Vector{<:AbstractFloat}[, unique_rxns=true, max_molecularity=2])
    rd = RxData(sd::SpeciesData, reacs::Vector{Vector{String}}, prods::Vector{Vector{String}}, 
                rsys::Vector{Dict{String, Any}}, psys::Vector{Dict{String, Any}}, 
                dH::Vector{<:AbstractFloat}, level::Int[, unique_rxns=true, max_molecularity=2])

Data container for reactions.

Constructor creates a new `RxData` reaction data store from 
lists of reactant and product SMILES by cross-referencing species
IDs with `sd`. `sd` should therefore already have all
species input here loaded in.

`reacs` and `prods` should be the raw SMILES
arrays from `ingest_cde_run()`, i.e. without any duplicate
species removed due to stoichiometry. This constructor will
determine stoichiometry and output the unique form of each
set of reactants/products.

Can optionally be given a `level` argument, indicating the 
exploration level which a set of reactions was discovered. 
If this is not provided, assumes reactions are entering the struct
at level 1.

By default, adds a maximum of 1 of each reaction type when
`unique_rxns = true`, and only admits reactions with a maximum
molecularity of 2 (i.e. bimolecular reactions).

Contains fields for:
* Number of reactions encountered (`nr`)
* Atom-mapped reaction SMILES for unambiguous linking of atom indices in reactants and products (`mapped_rxns`)
* Unique IDs of reactants for each reaction (`id_reacs`)
* Unique IDs of products for each reaction (`id_prods`)
* Stoichiometries of reactants for each reaction (`stoic_reacs`)
* Stoichiometries of products for each reaction (`stoic_prods`)
* Reaction enthalpies (`dH`)
* Reaction hashes, used for unique identification (`rhash`)
* The exploration level on which reactions were first found (`level_found`)
"""
function RxData(sd::SpeciesData{iType}, 
        reacs::Vector{Vector{String}}, prods::Vector{Vector{String}},
        rsys::Vector{Dict{String, Any}}, psys::Vector{Dict{String, Any}},
        dH::Vector{fType}, level::Int; 
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

    levels_final = [level for _ in 1:nr]
    return RxData{iType, fType}(nr, rxns_final, id_reacs_final, id_prods_final,
            stoic_reacs_final, stoic_prods_final, dH_final, hashes_final, levels_final)
end

function RxData(sd::SpeciesData{iType}, 
    reacs::Vector{Vector{String}}, prods::Vector{Vector{String}},
    rsys::Vector{Dict{String, Any}}, psys::Vector{Dict{String, Any}},
    dH::Vector{fType}; 
    unique_rxns=true, max_molecularity=2) where {iType, fType <: AbstractFloat}

    return RxData(sd, reacs, prods, rsys, psys, dH, 1; unique_rxns=unique_rxns, max_molecularity=max_molecularity)
end

"""
    push!(rd::RxData, sd::SpeciesData, reacs::Vector{Vector{String}}, prods::Vector{Vector{String}}, 
          rsys::Vector{Dict{String, Any}}, psys::Vector{Dict{String, Any}}, 
          dH::Vector{<:AbstractFloat}[, unique_rxns=true, max_molecularity=2])
    push!(rd::RxData, sd::SpeciesData, reacs::Vector{Vector{String}}, prods::Vector{Vector{String}}, 
          rsys::Vector{Dict{String, Any}}, psys::Vector{Dict{String, Any}}, 
          dH::Vector{<:AbstractFloat}, level::Int[, unique_rxns=true, max_molecularity=2])

Adds an array of reactions to `rd`.

Extends an `RxData` reaction data store from lists of reactant
and product SMILES by cross-referencing species IDs with `sd`.
`sd` should therefore already have all species input here loaded
in.

`reacs` and `prods` should be the raw SMILES
arrays from `ingest_cde_run()`, i.e. without any duplicate
species removed due to stoichiometry. This function will
determine stoichiometry and output the unique form of each
set of reactants/products.

Can optionally be given a `level` argument, indicating the 
exploration level which a set of reactions was discovered. 
If this is not provided, assumes reactions are entering the struct
at level 1.

By default, adds a maximum of 1 of each reaction type when
`unique_rxns = true`, and only admits reactions with a maximum
molecularity of 2 (i.e. bimolecular reactions).
"""
function Base.push!(rd::RxData{iType, fType}, sd::SpeciesData,
        reacs::Vector{Vector{String}}, prods::Vector{Vector{String}}, 
        rsys::Vector{Dict{String, Any}}, psys::Vector{Dict{String, Any}},
        dH::Vector{fType}, level::Int; 
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
            push!(rd.level_found, level)
            
        else
            dcounter += 1
        end
    end

    @debug " - $dcounter duplicate and $invcounter invalid reactions found."

    return
end

function Base.push!(rd::RxData{iType, fType}, sd::SpeciesData,
    reacs::Vector{Vector{String}}, prods::Vector{Vector{String}}, 
    rsys::Vector{Dict{String, Any}}, psys::Vector{Dict{String, Any}},
    dH::Vector{fType}; 
    unique_rxns=true, max_molecularity=2) where {iType, fType <: AbstractFloat}

    push!(rd, sd, reacs, prods, rsys, psys, dH, 1; unique_rxns=unique_rxns, max_molecularity=max_molecularity)
    return
end


"""
    get_reverse_rhash(sd, rd, rid)

Returns the reverse reaction hash for the reaction at `rid` in `rd`.

Useful when needing to identify if a reverse reaction
is already in a CRN without having to look through many
species permutations.
"""
function get_reverse_rhash(sd::SpeciesData, rd::RxData, rid)
    reacs = []
    for (i, sid) in enumerate(rd.id_reacs[rid])
        for _ in 1:rd.stoic_reacs[rid][i] 
            push!(reacs, sd.toStr[sid])
        end
    end
    sort!(reacs)
    prods = []
    for (i, sid) in enumerate(rd.id_prods[rid])
        for _ in 1:rd.stoic_prods[rid][i] 
            push!(prods, sd.toStr[sid])
        end
    end
    sort!(prods)

    forw_rhash = stable_hash(vcat(reacs, prods))
    @assert rd.rhash[rid] == forw_rhash
    rev_rhash = stable_hash(vcat(prods, reacs))
    return rev_rhash
end


"""
    init_network([iType=Int64, fType=Float64])

Initialises an empty reaction network.

Returns an empty `SpeciesData{iType}` and RxData{iType, fType}.
"""
function init_network(; iType=Int64, fType=Float64)
    sd = SpeciesData{iType}(
        Dict{String, iType}(), Dict{iType, String}(),
        0, 
        Dict{iType, Dict{String, Any}}(), 
        Dict{iType, Int}(), Dict()
    )
    rd = RxData{iType, fType}(
        0, String[],
        Vector{iType}[], Vector{iType}[], 
        Vector{iType}[], Vector{iType}[], 
        fType[], Vector{UInt8}[], Vector{Int}[]
    )

    return sd, rd
end


"""
    splice!(rd::RxData, rids::Vector{Int})

Removes reactions at indices `rids` from `rd`.
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