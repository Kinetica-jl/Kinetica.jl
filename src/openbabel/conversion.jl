"""
    smi_list, xyz_list = ingest_xyz_system(xyz_file[, fix_radicals])

Converts a molecule system from a single XYZ file into a list of SMILES strings and their respective ExtXYZ representataions.

OBCR can be used to attempt to fix Openbabel's radical structure
by enabling `fix_radicals`.
"""
function ingest_xyz_system(xyz_file::String; fix_radicals::Bool=true)
    pbmol = collect(pybel.readfile("xyz", xyz_file))[1]
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

    # Collate individual xyzs and convert to ExtXYZ frames.
    all_xyzs = ""
    for frag in fragments
        all_xyzs *= frag.write("xyz")
    end
    xyz_list = xyz_to_frames(all_xyzs)

    return smi_list, xyz_list
end


"""
    frames = xyz_to_frames(xyz)

Convenient handler for converting string-form xyz to ExtXYZ frames.

ExtXYZ is very resistant to reading in frames from memory,
instead requiring a file to load from. Handles creation of
a temporary file for this purpose, and reading back in of
the resulting frames.

Can also be used to go directly from a pbmol to frame:

    frames = xyz_to_frames(pbmol.write("xyz"))
"""
function xyz_to_frames(xyz::String)
    path, io = mktemp()
    write(io, xyz)
    close(io)
    return read_frames(path)
end