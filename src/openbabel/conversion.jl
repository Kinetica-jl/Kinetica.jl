"""
    ingest_xyz_system(xyz_str::String[, surfacefinder::Py, fix_radicals=true])

Converts a molecule system from a single XYZ string into a list of SMILES strings and their respective ExtXYZ representations.

If `surfaces` are passed, runs graph isomorphism checks on all
disconnected species for each surface's unit cell, isolating the
correct surface and absorption point to be used in the generated
species' SMILES.

OBCR can be used to attempt to fix Openbabel's radical structure
by enabling `fix_radicals`. If parsing directly from an XYZ file
is required, this can be achieved with

    smi_list, xyz_list = ingest_xyz_system(xyz_file_to_str(xyz_file))
"""
function ingest_xyz_system(xyz_str::String, surfdata::SurfaceData; fix_radicals=true)
    surf_elems = surfdata.finder.elements
    frame = xyz_to_frame(xyz_str)
    ads_slab = frame_to_atoms(frame)

    # If no surface is found, proceed with normal ingest.
    has_surface = pyconvert(bool, asesf.utils.has_elems(ads_slab, surf_elems))
    if !(has_surface)
        smi_list, xyz_list = ingest_xyz_system(xyz_str; fix_radicals)
        return smi_list, xyz_list
    end

    # Otherwise, separate adsorbates and predict surface sites.
    _, ads_molecules, sf_labels_per_molecule = surfdata.finder.predict(ads_slab)

    smi_list = []
    xyz_list = []
    for (ads_molecule, sf_labels) in zip(ads_molecules, sf_labels_per_molecule)
        ads_atom_idxs = pyconvert(Vector{Int}, sf_labels.keys())
        ads_frame = atoms_to_frame(ads_molecule)
        push!(xyz_list, ads_frame)
        ads_pbmol = pybel.readstring("xyz", frame_to_xyz(ads_frame))
        ads_smi = String(strip(pyconvert(String, ads_pbmol.write("can")), ['\n', '\t']))

        if fix_radicals && pyconvert(Bool, obcr.is_radical(ads_smi))
            ads_pbmol = obcr.fix_radicals(ads_pbmol)
            ads_pbmol.addh()
            ads_smi = String(strip(pyconvert(String, ads_pbmol.write("can")), ['\n', '\t']))
        end

        elem_replacements = []
        site_atomic_number = 100
        for atom_idx in ads_atom_idxs
            # Determine SMILES label for surface site that adsorbed atom is on.
            label = pyconvert(String, sf_labels[atom_idx]["site"])
            split_labels = split(label, '_')
            surf_label = join(split_labels[begin:end-1], '_')
            site_label = String(split_labels[end])
            surf_idx = surfdata.nameToInt[surf_label]
            site_idx = surfdata.surfaces[surf_idx].siteids[site_label]
            smi_label = "X$(surf_idx)_$(site_idx)"

            # Attach unique atoms of atomic number 100+ in OB to adsorbed atom
            site_elem = pyconvert(String, pybel.ob.GetSymbol(site_atomic_number))
            site_atom = ads_pbmol.OBMol.NewAtom()
            site_atom.SetAtomicNum(site_atomic_number)
            ads_pbmol.OBMol.AddBond(atom_idx+1, site_atom.GetIdx(), 1)
            
            push!(elem_replacements, Pair(site_elem, smi_label))
        end

        # Generate new SMILES, replace dummy atoms with ads labels.
        ads_smi = String(strip(pyconvert(String, ads_pbmol.write("can")), ['\n', '\t']))
        replace(ads_smi, elem_replacements...)

        push!(smi_list, ads_smi)
    end

    return smi_list, xyz_list
end


function ingest_xyz_system(xyz_str::String; fix_radicals=true)
    pbmol = pybel.readstring("xyz", xyz_str)
    fragments = [pybel.Molecule(obmol) for obmol in pbmol.OBMol.Separate()]
    n = length(fragments)
    smi_list = String[String(strip(pyconvert(String, frag.write("can")), ['\n', '\t'])) for frag in fragments]

    # Fix radical structures if requested.
    if fix_radicals
        for i in 1:n
            if pyconvert(Bool, obcr.is_radical(smi_list[i]))
                fragments[i] = obcr.fix_radicals(fragments[i])
                fragments[i].addh()
                smi_list[i] = String(strip(pyconvert(String, fragments[i].write("can")), ['\n', '\t']))
            end
        end
    end

    # Convert to ExtXYZ frames.
    xyz_list = [xyz_to_frame(pyconvert(String, frag.write("xyz"))) for frag in fragments]

    return smi_list, xyz_list
end


"""
    xyz_to_frame(xyz::String)

Convenient handler for converting string-form xyz to ExtXYZ frame.

Reformats string-form xyz into an `IOBuffer` so it can
be read in-memory by `read_frames`.

Can also be used to go directly from a Pybel pbmol to frame:

    frame = xyz_to_frame(pbmol.write("xyz"))
"""
function xyz_to_frame(xyz::String)
    iob = IOBuffer(xyz)
    frame = read_frame(iob)
    return frame
end


function _py_xyz_to_frame(xyz::String)
    function map_dtypes(descr)
        out = []
        for dtype in descr
            pytype = pyconvert(String, dtype[1])
            jltype = nothing
            if 'U' in pytype
                jltype = String
            elseif 'f' in pytype
                jltype = Float64
            elseif 'i' in pytype
                jltype = Int
            else
                error("Unknown dtype in ExtXYZ arrays: $(pytype)")
            end
            if pylen(dtype) == 3
                jltype = Vector{jltype}
            end
            push!(out, jltype)
        end
        return out
    end

    iob = IOBuffer(xyz)
    na, info, arrays, _ = pyextxyz.extxyz.read_frame_dicts(split(String(take!(iob)), '\n'), use_regex=false)
    close(iob)

    # Arrays in info dict need to be pure Julia for ExtXYZ.jl writes,
    # but also needs to be a Dict{String, Any}.
    info = pyconvert(Dict{String, Array}, info)
    info = Dict{String, Any}(k => v for (k, v) in pairs(info))
    
    na = pyconvert(Int, na)
    array_keys = pyconvert(Vector{String}, arrays.dtype.names)
    dtype_map = map_dtypes(arrays.dtype.descr)

    jlarrays = Dict{String, Any}()
    for (i, (key, dtype)) in enumerate(zip(array_keys, dtype_map))
        if dtype <: Vector
            if na == 1
                jlarrays[key] = reduce(hcat, [pyconvert(dtype, arrays.item()[i-1])]) 
            else
                jlarrays[key] = reduce(hcat, [pyconvert(dtype, a[i-1]) for a in arrays])
            end
        else
            if na == 1
                jlarrays[key] = [pyconvert(dtype, arrays.item()[i-1])]
            else
                jlarrays[key] = [pyconvert(dtype, a[i-1]) for a in arrays]
            end
        end
    end

    frame = Dict{String, Any}("N_atoms" => na, "info" => info, "arrays" => jlarrays)
    return frame
end


"""
    frame_to_xyz(frame::Dict{String, Any})

Converts an ExtXYZ frame into string-form XYZ.

Does not convert frames to XYZ files. If this is required, use
the `write_frames` function of ExtXYZ.
"""
function frame_to_xyz(frame::Dict{String, Any})
    na = frame["N_atoms"]
    s = "$na\n"
    comment = join(["$key=$value" for (key, value) in frame["info"]], " ")
    s *= "$comment\n"
    for i in 1:na
        s *= "$(frame["arrays"]["species"][i]) $(frame["arrays"]["pos"][1, i]) $(frame["arrays"]["pos"][2, i]) $(frame["arrays"]["pos"][3, i])\n"
    end
    return s
end


"""
    xyz_from_smiles(smi::String[, generator::Symbol=openbabel, seed=-1])
    xyz_from_smiles(smi::String, saveto::String[, generator::Symbol=:openbabel,  overwrite=true, seed=-1])

Generates approximate coordinates for a molecule from SMILES.

If provided, saves to an XYZ file at `saveto`. If already present,
this file will only be overwritten when `overwrite=true` Otherwise
returns the string-form XYZ without saving to file.

Can generate geometries using either Openbabel or RDKit, chosen
using `generator=:openbabel` or `generator=:rdkit`. RDKit-based
generation can be seeded for consistency using the `seed` argument,
while trying to supply a seed to Openbabel will throw an error.

Does not currently support the creation of surface-bound species from
SMILES, although this feature is planned.
"""
xyz_from_smiles(smi::String; generator::Symbol=:openbabel, seed=-1) = xyz_from_smiles(Val(generator), SpeciesStyle(smi), smi, seed) 
function xyz_from_smiles(::Val{:openbabel}, ::GasSpecies, smi::String, seed)
    if seed != -1
        throw(ArgumentError("OpenBabel generator does not support seeded generation of molecular geometries."))
    end
    pbmol = pybel.readstring("smi", smi)
    ff = smi in ["[H][H]", "[H]"] ? "uff" : "mmff94"
    pbmol.make3D(forcefield=ff)
    return pyconvert(String, pbmol.write("xyz"))
end
function xyz_from_smiles(::Val{:openbabel}, ::SurfaceSpecies, smi::String, seed)
    throw(ArgumentError("Kinetica does not currently support creation of surface-bound species from SMILES."))
end

xyz_from_smiles(smi::String, saveto::String; generator::Symbol=:openbabel, overwrite::Bool=true, seed=-1) = xyz_from_smiles(Val(generator), SpeciesStyle(smi), smi, saveto, overwrite, seed)
function xyz_from_smiles(::Val{:openbabel}, ::GasSpecies, smi::String, saveto::String, overwrite::Bool, seed)
    if seed != -1
        throw(ArgumentError("OpenBabel generator does not support seeded generation of molecular geometries."))
    end
    pbmol = pybel.readstring("smi", smi)
    ff = smi in ["[H][H]", "[H]"] ? "uff" : "mmff94"
    pbmol.make3D(forcefield=ff)
    pbmol.write("xyz", saveto; overwrite=overwrite)
    return
end
function xyz_from_smiles(::Val{:openbabel}, ::SurfaceSpecies, smi::String, saveto::String, overwrite::Bool, seed)
    throw(ArgumentError("Kinetica does not currently support creation of surface-bound species from SMILES."))
end


"""
    frame_from_smiles(smi::String)

Generates approximate coordinates for a molecule from SMILES.

Uses `xyz_from_smiles` through the Openbabel generator
internally. If greater control over generation is
required, call this separately. Returns an ExtXYZ
frame dictionary.
"""
function frame_from_smiles(smi::String)
    frame = xyz_to_frame(xyz_from_smiles(smi))
    return frame
end


"""
    xyz_file_to_str(xyz_file::String)

Converts contents of an XYZ file into string-form XYZ.
"""
function xyz_file_to_str(xyz_file::String)
    return read(xyz_file, String)
end