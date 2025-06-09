function xyz_from_smiles(::Val{:rdkit}, ::GasSpecies, smi::String, seed)
    mol = rdChem.MolFromSmiles(smi)
    mol = rdChem.AddHs(mol)
    rdChem.rdDistGeom.EmbedMolecule(mol, randomSeed=seed)
    if smi in ["[H][H]", "[H]"]
        rdChem.rdForceFieldHelpers.MMFFOptimizeMolecule(mol)
    else
        rdChem.rdForceFieldHelpers.UFFOptimizeMolecule(mol)
    end

    return pyconvert(String, rdChem.MolToXYZBlock(mol))
end
function xyz_from_smiles(::Val{:rdkit}, ::SurfaceSpecies, smi::String, seed)
    throw(ArgumentError("Kinetica does not currently support creation of surface-bound species from SMILES."))
end

function xyz_from_smiles(::Val{:rdkit}, ::GasSpecies, smi::String, saveto::String, overwrite::Bool, seed)
    mol = rdChem.MolFromSmiles(smi)
    mol = rdChem.AddHs(mol)
    rdChem.rdDistGeom.EmbedMolecule(mol, randomSeed=seed)
    if smi in ["[H][H]", "[H]"]
        rdChem.rdForceFieldHelpers.MMFFOptimizeMolecule(mol)
    else
        rdChem.rdForceFieldHelpers.UFFOptimizeMolecule(mol)
    end

    if !overwrite throw(ArgumentError(
        "RDKit generator always overwrites previous files, please save XYZ to a different location."
    )) end
    rdChem.MolToXYZFile(mol, saveto)
    return
end
function xyz_from_smiles(::Val{:rdkit}, ::SurfaceSpecies, smi::String, saveto::String, overwrite::Bool, seed)
    throw(ArgumentError("Kinetica does not currently support creation of surface-bound species from SMILES."))
end

"""
    frame_to_rdkit(frame::Dict{String, Any}[, with_coords=false])

Converts an ExtXYZ frame to an RDKit `Mol` object.

Since Rdkit often fails to percieve bonding from geometry alone,
this creates an explicitly single-bonded system with no perception
of bond orders or unpaired electrons. It therefore functions as
a raw connectivity map.

If bond order determination is required, atoms can be given their
3D coordinates using `with_coords`, such that bond orders can be
inferred by bond lengths.
"""
function frame_to_rdkit(frame::Dict{String, Any}; with_coords=false)
    pbmol = pybel.readstring("xyz", frame_to_xyz(frame))
    rdmol = rdChem.rdchem.EditableMol(rdChem.rdchem.Mol())
    for i in 1:frame["N_atoms"]
        rdatom = rdChem.rdchem.Atom(frame["arrays"]["species"][i])
        rdatom.SetAtomMapNum(i)
        rdmol.AddAtom(rdatom)
    end

    rdmol = frame_to_rdkit_remap_atoms(pbmol, rdmol)

    if with_coords
        conf = rdChem.Conformer(frame["N_atoms"])
        rdmol.AddConformer(conf)
        conf = rdmol.GetConformer()
        for i in 1:frame["N_atoms"]
            conf.SetAtomPosition(i-1, rdGeometry.Point3D(frame["arrays"]["pos"][:, i]...))
        end
    end

    return rdmol
end


"""
    atom_map_smiles(frame::Dict{String, Any}, smi::String[, allow_mismatch=false])
    atom_map_smiles(frame::Dict{String, Any}, smi::String[, resub_ads_bonding=false])

Maps the atom indices from `frame` to the atoms in `smi`.

Uses a raw connectivity map (purely single-bonded) representation
of the molecule system represented in both the ExtXYZ geometry `frame`
and its canonical SMILES `smi` to construct a substructure match, 
which can be used to map atom indices to every atom in `smi`.

Useful when atom-mapped geometries of a reaction's reactants and
products are accessible, as by calling this function for each, a
fully atom-mapped reaction SMILES can be constructed, allowing for
later reconstruction of atom-mapped geometries. Permits imperfect
matches between representations if `allow_mismatch` is set to
`true` (not recommended).

When given a surface SMILES and a corresponding adsorbate frame
without any surface atoms, makes a set of substitutions that
allows for atom mapping the adsorbate without assigning any
indices to the surface sites themselves. In this case, the
`allow_ads_mismatch` argument is disabled, as there must always
be a slight mismatch to account for these substitutions. This
method allows for the coordination of sites to be passed to the
atom-mapped SMILES with `resub_ads_bonding=true`.

Heavily based on the implementation in Colin Grambow's `ard_gsm`
package: https://github.com/cgrambow/ard_gsm/tree/v1.0.0
"""
atom_map_smiles(frame::Dict{String, Any}, smi::String; kwargs...) = atom_map_smiles(XYZStyle(frame), SpeciesStyle(smi), frame, smi; kwargs...)
function atom_map_smiles(::FreeXYZ, ::GasSpecies, frame::Dict{String, Any}, smi::String; allow_mismatch=false)
    atoms_in_mol_true = Dict{String, Int}()
    for i in 1:frame["N_atoms"]
        elem = frame["arrays"]["species"][i]
        atoms_in_mol_true[elem] = get(atoms_in_mol_true, elem, 0) + 1
    end

    mol_sanitised = rdChem.MolFromSmiles(smi)
    if pyconvert(Bool, mol_sanitised == pybuiltins.None)
        throw(ErrorException("Unable to create an RDKit mol from provided SMILES."))
    end
    mol_sanitised = rdChem.AddHs(mol_sanitised)
    atoms_in_mol_sani = Dict{String, Int}()
    for atom in mol_sanitised.GetAtoms()
        elem = pyconvert(String, atom.GetSymbol())
        atoms_in_mol_sani[elem] = get(atoms_in_mol_sani, elem, 0) + 1
    end
    if atoms_in_mol_true != atoms_in_mol_sani && !(allow_mismatch)
        println(smi)
        println(frame)
        println(atoms_in_mol_true)
        println(atoms_in_mol_sani)
        throw(ErrorException("Unable to match SMILES atoms to XYZ atoms."))
    end

    mol_with_map = frame_to_rdkit(frame)
    mol_sani_sb = rdChem.Mol(mol_sanitised)
    for bond in mol_sani_sb.GetBonds()
        bond.SetBondType(rdChem.rdchem.BondType."SINGLE")
    end

    match = pyconvert(Vector, mol_sani_sb.GetSubstructMatch(mol_with_map))
    if pyconvert(Int, mol_with_map.GetNumAtoms()) != length(match)
        println(mol_with_map.GetNumAtoms())
        throw(ErrorException("Incorrect number of atoms when matching substruct during atom mapping."))
    end
    for atom in mol_with_map.GetAtoms()
        idx = match[pyconvert(Int, atom.GetIdx()) + 1]
        map_num = atom.GetAtomMapNum()
        mol_sanitised.GetAtomWithIdx(idx).SetAtomMapNum(map_num)
    end

    # Disable any dative bonds that RDKit wants to make, since
    # OpenBabel can't handle them.
    for bond in mol_sanitised.GetBonds()
        if pyis(bond.GetBondType(), rdChem.rdchem.BondType."DATIVE")
            bond.SetBondType(rdChem.rdchem.BondType."SINGLE")
        end
    end

    return pyconvert(String, rdChem.MolToSmiles(mol_sanitised))
end
function atom_map_smiles(::Union{AdsorbateXYZ, OnSurfaceXYZ}, ::GasSpecies, ::Dict{String, Any}, ::String; kwargs...)
    throw(ErrorException("Unable to map gas-phase SMILES from a surface-bound geometry."))
end

function atom_map_smiles(::AdsorbateXYZ, ::SurfaceSpecies, frame::Dict{String, Any}, smi::String; resub_ads_bonding=false)
    surfid = get_surfid(smi)
    siteids = get_surf_siteids(smi)
    pt = rdChem.GetPeriodicTable()

    # Substitute out bonding to surface sites.
    site_replacements = []
    re = r"(?<=[#=])(\[X\d_\d\])" # Search for #/= before site tags
    m = match(re, smi)
    while !isnothing(m)
        push!(site_replacements, Pair(smi[m.offset-1]*m.match, m.match))
        m = match(re, smi, m.offset+length(m.match))
    end
    re = r"(\[X\d_\d\])(?=[#=])" # Search for #/= after site tags
    m = match(re, smi)
    while !isnothing(m)
        push!(site_replacements, Pair(m.match*smi[m.offset+length(m.match)], m.match))
        m = match(re, smi, m.offset+length(m.match)+1)
    end
    smi_subbed = replace(smi, site_replacements...)

    # Substitute surface site tags with unique elements.
    site_atomic_number = 100
    elem_replacements = []
    for siteid in siteids
        elem = pyconvert(String, pt.GetElementSymbol(site_atomic_number))
        push!(elem_replacements, Pair("X$(surfid)_$(siteid)", elem))
        site_atomic_number += 1
    end
    smi_replaced = replace(smi_subbed, elem_replacements...)

    # Generate atom mapped SMILES, ignoring substituted elements.
    amsmi = atom_map_smiles(FreeXYZ(), GasSpecies(), frame, smi_replaced; allow_mismatch=true)

    # Substitute site tags back into atom mapped SMILES.
    elem_replacements_flipped = [Pair(e[2], e[1]) for e in elem_replacements]
    amsmi_replaced = replace(amsmi, elem_replacements_flipped...)

    # Substitute bonding back onto surface sites.
    if resub_ads_bonding
        site_replacements_flipped = [Pair(e[2], e[1]) for e in site_replacements]
        amsmi_final = replace(amsmi_replaced, site_replacements_flipped...)
    else
        amsmi_final = amsmi_replaced
    end

    return amsmi_final
end
function atom_map_smiles(::FreeXYZ, ::SurfaceSpecies, frame::Dict{String, Any}, smi::String; resub_ads_bonding=false)
    @warn "Attempting to map an unadsorbed geometry to a surface SMILES. This is not supported and may fail."
    return atom_map_smiles(AdsorbateXYZ(), SurfaceSpecies(), frame, smi; resub_ads_bonding=resub_ads_bonding)
end
function atom_map_smiles(::OnSurfaceXYZ, ::SurfaceSpecies, ::Dict{String, Any}, ::String; kwargs...)
    throw(ErrorException("Unable to generate mapped surface SMILES from a surface-bound species. Remove surface atoms from geometry and pass isolated adsorbate."))
end


"""
    atom_map_frame(am_smi::String, frame::Dict{String, Any})

Maps the atom indices from `am_smi` to the geometry in `frame`.

Uses an atom-mapped SMILES `am_smi` to reorder atoms in an
ExtXYZ system geometry `frame` by substructure matching, as
in `atom_map_smiles`.

Made a little more complicated than generating an atom-mapped 
SMILES from a geometry by the fact the Rdkit can't seem to output
XYZs using internal atom maps, and correctly binding new custom
conformers to existing `Mol`s without calling `EmbedMolecule` is
not simple. Instead, constructs a transfer array for atom indices
and applies this directly to the `frame`.
"""
atom_map_frame(am_smi::String, frame::Dict{String, Any}) = atom_map_frame(SpeciesStyle(am_smi), XYZStyle(frame), am_smi, frame)
function atom_map_frame(::GasSpecies, ::FreeXYZ, am_smi::String, frame::Dict{String, Any}; allow_mismatch=false, ignore_elements=nothing)
    smiles_params = rdChem.SmilesParserParams()
    smiles_params.removeHs = false
    smiles_params.sanitize = false
    mol_template = rdChem.MolFromSmiles(am_smi, smiles_params)
    atoms_in_mol_template = Dict{String, Int}()
    for atom in mol_template.GetAtoms()
        elem = pyconvert(String, atom.GetSymbol())
        atoms_in_mol_template[elem] = get(atoms_in_mol_template, elem, 0) + 1
    end
    for bond in mol_template.GetBonds()
        bond.SetBondType(rdChem.rdchem.BondType."SINGLE")
    end

    mol_target = frame_to_rdkit(frame)
    atoms_in_mol_target = Dict{String, Int}()
    for atom in mol_target.GetAtoms()
        elem = pyconvert(String, atom.GetSymbol())
        atoms_in_mol_target[elem] = get(atoms_in_mol_target, elem, 0) + 1
    end
    if atoms_in_mol_template != atoms_in_mol_target && !(allow_mismatch)
        println(amsmi)
        println(frame)
        println(atoms_in_mol_template)
        println(atoms_in_mol_target)
        throw(ErrorException("Unable to match SMILES atoms to XYZ atoms."))
    end

    mol_target_sb = rdChem.Mol(mol_target)
    for bond in mol_target_sb.GetBonds()
        bond.SetBondType(rdChem.rdchem.BondType."SINGLE")
    end
    for atom in mol_target_sb.GetAtoms()
        atom.SetAtomMapNum(0)
    end

    if !isnothing(ignore_elements)
        remove_idxs = []
        for atom in mol_template.GetAtoms()
            elem = pyconvert(String, atom.GetSymbol())
            if elem in ignore_elements
                push!(remove_idxs, pyconvert(Int, atom.GetIdx()))
            end
        end
        if length(remove_idxs) > 0
            mol_template = rdChem.RWMol(mol_template)
            offset = 0
            for idx in remove_idxs
                mol_template.RemoveAtom(idx-offset)
                offset += 1
            end
        end
    end
    
    match = pyconvert(Vector, mol_target_sb.GetSubstructMatch(mol_template))
    if pyconvert(Int, mol_template.GetNumAtoms()) != length(match)
        throw(ErrorException("Incorrect number of atoms when matching substruct during atom mapping ($(length(match)) matched, $(mol_template.GetNumAtoms()) expected)."))
    end
    for atom in mol_template.GetAtoms()
        idx = match[pyconvert(Int, atom.GetIdx()) + 1]
        map_num = atom.GetAtomMapNum()
        mol_target.GetAtomWithIdx(idx).SetAtomMapNum(map_num)
    end
    transfer = zeros(Int, frame["N_atoms"])
    for atom in mol_target.GetAtoms()
        transfer[pyconvert(Int, atom.GetIdx())+1] = pyconvert(Int, atom.GetAtomMapNum())
    end

    new_frame = deepcopy(frame)
    for (i, m) in enumerate(transfer)
        new_frame["arrays"]["pos"][:, m] = frame["arrays"]["pos"][:, i]
        new_frame["arrays"]["species"][m] = frame["arrays"]["species"][i]
    end

    return new_frame
end
function atom_map_frame(::GasSpecies, ::Union{AdsorbateXYZ, OnSurfaceXYZ}, am_smi::String, frame::Dict{String, Any})
    throw(ErrorException("Unable to map a surface-bound geometry from a gas-phase SMILES."))
end

function atom_map_frame(::SurfaceSpecies, ::AdsorbateXYZ, am_smi::String, frame::Dict{String, Any})
    surfid = get_surfid(am_smi)
    siteids = get_surf_siteids(am_smi)
    pt = rdChem.GetPeriodicTable()

    # Substitute out bonding to surface sites.
    site_replacements = []
    re = r"(?<=[#=])(\[X\d_\d\])" # Search for #/= before site tags
    m = match(re, am_smi)
    while !isnothing(m)
        push!(site_replacements, Pair(am_smi[m.offset-1]*m.match, m.match))
        m = match(re, am_smi, m.offset+length(m.match))
    end
    re = r"(\[X\d_\d\])(?=[#=])" # Search for #/= after site tags
    m = match(re, am_smi)
    while !isnothing(m)
        push!(site_replacements, Pair(m.match*am_smi[m.offset+length(m.match)], m.match))
        m = match(re, am_smi, m.offset+length(m.match)+1)
    end
    amsmi_subbed = replace(am_smi, site_replacements...)

    # Substitute surface site tags with unique elements.
    site_atomic_number = 100
    elem_replacements = []
    remove_elems = []
    for siteid in siteids
        elem = pyconvert(String, pt.GetElementSymbol(site_atomic_number))
        push!(elem_replacements, Pair("X$(surfid)_$(siteid)", elem))
        push!(remove_elems, elem)
        site_atomic_number += 1
    end
    amsmi_replaced = replace(amsmi_subbed, elem_replacements...)

    # Atom map frame from SMILES, ignoring substituted elements.
    amframe = atom_map_frame(GasSpecies(), FreeXYZ(), amsmi_replaced, frame; 
                             allow_mismatch=true, ignore_elements=remove_elems)

    return amframe
end
function atom_map_frame(::SurfaceSpecies, ::OnSurfaceXYZ, am_smi::String, frame::Dict{String, Any})
    if !haskey(frame["arrays"], "tags")
        throw(ErrorException("Unable to map surface-bound geometry - missing surface-adsorbate tags."))
    end
    # Remove surface atoms from geometry
    ads_idxs = findall(frame["arrays"]["tags"] .== 0)
    ads_frame = Dict{String, Any}("arrays" => Dict{String, Any}(), "info" => Dict{String, Any}())
    ads_frame["N_atoms"] = length(ads_idxs)
    ads_frame["arrays"]["pos"] = frame["arrays"]["pos"][:, ads_idxs]
    ads_frame["arrays"]["species"] = frame["arrays"]["species"][ads_idxs]
    ads_frame["info"]["adsorbate"] = "true"

    # Separately map adsorbate to SMILES.
    ads_frame_mapped = atom_map_frame(SurfaceSpecies(), AdsorbateXYZ(), am_smi, ads_frame)

    # Recombine mapped adsorbate and surface.
    mapped_frame = deepcopy(frame)
    surf_idxs = findall(frame["arrays"]["tags"] .!= 0)
    surf_tags = frame["arrays"]["tags"][surf_idxs]
    surf_pos = frame["arrays"]["pos"][:, surf_idxs]
    surf_species = frame["arrays"]["species"][surf_idxs]
    mapped_frame["arrays"]["pos"] = hcat(ads_frame_mapped["arrays"]["pos"], surf_pos)
    mapped_frame["arrays"]["species"] = vcat(ads_frame_mapped["arrays"]["species"], surf_species)
    mapped_frame["arrays"]["tags"] = vcat(zeros(Int, length(ads_idxs)), surf_tags)

    return mapped_frame
end
function atom_map_frame(::SurfaceSpecies, ::FreeXYZ, am_smi::String, frame::Dict{String, Any})
    throw(ErrorException("Unable to map a gas-phase geometry from a surface-bound SMILES."))
end