"""
    ingest_xyz_system(xyz_str::String[, fix_radicals=true])

Converts a molecule system from a single XYZ string into a list of SMILES strings and their respective ExtXYZ representataions.

OBCR can be used to attempt to fix Openbabel's radical structure
by enabling `fix_radicals`. If parsing directly from an XYZ file
is required, this can be achieved with

    smi_list, xyz_list = ingest_xyz_system(xyz_file_to_str(xyz_file))
"""
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

ExtXYZ is very resistant to reading in frames from memory,
instead requiring a file to load from. Circumvents this by
falling back to a Python implementation of the reader which
can accept a list of lines from the file to iterate over.
Constructs this list within an in-memory IOBuffer to avoid
large IO overhead of writing to disk (repeatedly).

Can also be used to go directly from a Pybel pbmol to frame:

    frame = xyz_to_frame(pbmol.write("xyz"))
"""
function xyz_to_frame(xyz::String)
    iob = IOBuffer(xyz)
    na, info, arrays, _ = pyextxyz.extxyz.read_frame_dicts(split(String(take!(iob)), '\n'), use_regex=false)
    close(iob)

    info = pyconvert(Dict{String, Any}, info)
    na = pyconvert(Int, na)
    if na == 1
        arrays = Dict{String, Any}("species" => [pyconvert(String, arrays.item()[0])], "pos" => reduce(hcat, [pyconvert(Vector, arrays.item()[1])]))
    else
        arrays = Dict{String, Any}("species" => [pyconvert(String, a[0]) for a in arrays], "pos" => reduce(hcat, [pyconvert(Vector, a[1]) for a in arrays]))
    end
    frame = Dict{String, Any}("N_atoms" => na, "info" => info, "arrays" => arrays)
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
"""
xyz_from_smiles(smi::String; generator::Symbol=:openbabel, seed=-1) = xyz_from_smiles(Val(generator), smi, seed) 
function xyz_from_smiles(::Val{:openbabel}, smi::String, seed)
    if seed != -1
        throw(ArgumentError("OpenBabel generator does not support seeded generation of molecular geometries."))
    end
    pbmol = pybel.readstring("smi", smi)
    ff = smi in ["[H][H]", "[H]"] ? "uff" : "mmff94"
    pbmol.make3D(forcefield=ff)
    return pyconvert(String, pbmol.write("xyz"))
end

xyz_from_smiles(smi::String, saveto::String; generator::Symbol=:openbabel, overwrite::Bool=true, seed=-1) = xyz_from_smiles(Val(generator), smi, saveto, overwrite, seed)
function xyz_from_smiles(::Val{:openbabel}, smi::String, saveto::String, overwrite::Bool, seed)
    if seed != -1
        throw(ArgumentError("OpenBabel generator does not support seeded generation of molecular geometries."))
    end
    pbmol = pybel.readstring("smi", smi)
    ff = smi in ["[H][H]", "[H]"] ? "uff" : "mmff94"
    pbmol.make3D(forcefield=ff)
    pbmol.write("xyz", saveto; overwrite=overwrite)
    return
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
    pbmol = collect(pybel.readfile("xyz", xyz_file))[1]
    if pyconvert(String, pbmol._gettitle()) == xyz_file
        pbmol._settitle("")
    end
    return pyconvert(String, pbmol.write("xyz"))
end