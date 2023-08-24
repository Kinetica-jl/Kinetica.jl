"""
Bidirectional String-Int dictionary for chemical species.

Contains fields for:
* SMILES string -> integer ID dictionary (`toInt`)
* Integer ID -> SMILES string dictionary (`toStr`)
* Number of species (`n`)
* ExtXYZ structures of species (`xyz`)
* Dictionary of per-species cached values (`cache`)
"""
mutable struct SpeciesData
    toInt::Dict{String, Int}
    toStr::Dict{Int, String}
    n::Int
    xyz::Dict{Int, Dict{String, Any}}
    cache::Dict{Symbol, Any}
end

"""
    sd = SpeciesData(smi_list, xyz_list)

Outer constructor method for `SpeciesData`, allowing for construction
from a list of SMILES strings and ExtXYZ structures.
"""
function SpeciesData(smi_list, xyz_list)
    n = length(smi_list)
    if n == 0
        return SpeciesData(Dict(), Dict(), 0, Dict(), Dict())
    end

    return SpeciesData(
        Dict(smi => i for (i, smi) in enumerate(smi_list)),
        Dict(i => smi for (i, smi) in enumerate(smi_list)),
        n,
        Dict(i => x for (i, x) in enumerate(xyz_list)),
        Dict()
    )
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
function Base.push!(sd::SpeciesData, smi, xyz)
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
    push_unique!(sd, smi, xyz)

Add a species SMILES to a `SpeciesData`, as long as it does not already exist there.
"""
function push_unique!(sd::SpeciesData, smi, xyz)
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
