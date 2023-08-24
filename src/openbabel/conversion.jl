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
    smi_list = [strip(frag.write("can"), ['\n', '\t']) for frag in fragments]

    # Fix radical structures if requested.
    if fix_radicals
        for i in 1:n
            if obcr.is_radical(smi_list[i])
                fragments[i] = obcr.fix_radicals(fragments[i])
                fragments[i].addh()
                smi_list[i] = strip(fragments[i].write("can"), ['\n', '\t'])
            end
        end
    end

    # Pretty horrible way of converting to ExtXYZ dicts, but read_frames 
    # only accepts actual files it seems.
    all_xyzs = ""
    for frag in fragments
        all_xyzs *= frag.write("xyz")
    end
    open("/tmp/kinetica.xyz", "w") do f
        write(f, all_xyzs)
    end
    xyz_list = read_frames("/tmp/kinetica.xyz")
    rm("/tmp/kinetica.xyz")

    return smi_list, xyz_list
end