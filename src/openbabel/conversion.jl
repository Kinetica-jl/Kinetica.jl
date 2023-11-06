"""
    smi_list, xyz_list = ingest_xyz_system(xyz_str[, fix_radicals])

Converts a molecule system from a single XYZ string into a list of SMILES strings and their respective ExtXYZ representataions.

OBCR can be used to attempt to fix Openbabel's radical structure
by enabling `fix_radicals`. If parsing directly from an XYZ file
is required, this can be achieved with

    smi_list, xyz_list = ingest_xyz_system(xyz_file_to_str(xyz_file))
"""
function ingest_xyz_system(xyz_str::String; fix_radicals::Bool=true)
    pbmol = pybel.readstring("xyz", xyz_str)
    fragments = [pybel.Molecule(obmol) for obmol in pbmol.OBMol.Separate()]
    n = length(fragments)
    smi_list = String[String(strip(frag.write("can"), ['\n', '\t'])) for frag in fragments]

    # Fix radical structures if requested.
    if fix_radicals
        for i in 1:n
            if obcr.is_radical(smi_list[i])
                fragments[i] = obcr.fix_radicals(fragments[i])
                fragments[i].addh()
                smi_list[i] = String(strip(fragments[i].write("can"), ['\n', '\t']))
            end
        end
    end

    # Convert to ExtXYZ frames.
    xyz_list = [xyz_to_frame(frag.write("xyz")) for frag in fragments]

    return smi_list, xyz_list
end


"""
    frame = xyz_to_frame(xyz)

Convenient handler for converting string-form xyz to ExtXYZ frame.

ExtXYZ is very resistant to reading in frames from memory,
instead requiring a file to load from. Circumvents this by
falling back to a Python implementation of the reader which
can accept a list of lines from the file to iterate over.
Constructs this list within an in-memory IOBuffer to avoid
large IO overhead of writing to disk (repeatedly).

Can also be used to go directly from a pbmol to frame:

    frame = xyz_to_frame(pbmol.write("xyz"))
"""
function xyz_to_frame(xyz::String)
    iob = IOBuffer(xyz)
    na, info, arrays, _ = pyextxyz.extxyz.read_frame_dicts(split(String(take!(iob)), '\n'), use_regex=false)
    close(iob)

    info = Dict{String, Any}(key => val for (key, val) in info)
    if na == 1
        arrays = Dict{String, Any}("species" => [arrays.item()[1]], "pos" => reduce(hcat, [arrays.item()[2]]))
    else
        arrays = Dict{String, Any}("species" => [a[1] for a in arrays], "pos" => reduce(hcat, [a[2] for a in arrays]))
    end
    frame = Dict{String, Any}("N_atoms" => na, "info" => info, "arrays" => arrays)
    return frame
end


"""
    xyz_str = frame_to_xyz(frame)

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
    xyz_from_smiles(smi[, saveto, overwrite])

Generates approximate coordinates for a molecule from SMILES.

If provided, saves to an XYZ file at `saveto`. Otherwise returns
the string-form XYZ without saving to file.
"""
function xyz_from_smiles(smi::String)
    pbmol = pybel.readstring("smi", smi)
    ff = smi == "[H][H]" ? "uff" : "mmff94"
    pbmol.make3D(forcefield=ff)
    return pbmol.write("xyz")
end

function xyz_from_smiles(smi::String, saveto::String; overwrite::Bool=true)
    pbmol = pybel.readstring("smi", smi)
    ff = smi == "[H][H]" ? "uff" : "mmff94"
    pbmol.make3D(forcefield=ff)
    pbmol.write("xyz", saveto; overwrite=overwrite)
end


"""
    frame_from_smiles(smi)

Generates approximate coordinates for a molecule from SMILES.

Returns an ExtXYZ frame dictionary.
"""
function frame_from_smiles(smi::String)
    frame = xyz_to_frame(xyz_from_smiles(smi))
    return frame
end


"""
    xyz_str = xyz_file_to_str(xyz_file)

Converts contents of an XYZ file into string-form XYZ.
"""
function xyz_file_to_str(xyz_file::String)
    pbmol = collect(pybel.readfile("xyz", xyz_file))[1]
    return pbmol.write("xyz")
end