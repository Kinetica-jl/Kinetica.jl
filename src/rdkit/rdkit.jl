function xyz_from_smiles(::Val{:rdkit}, smi::String, seed)
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

function xyz_from_smiles(::Val{:rdkit}, smi::String, saveto::String, overwrite::Bool, seed)
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
    atom_map_smiles(frame::Dict{String, Any}, smi::String)

Maps the atom indices from `frame` to the atoms in `smi`.

Uses a raw connectivity map (purely single-bonded) representation
of the molecule system represented in both the ExtXYZ geometry `frame`
and its canonical SMILES `smi` to construct a substructure match, 
which can be used to map atom indices to every atom in `smi`.

Useful when atom-mapped geometries of a reaction's reactants and
products are accessible, as by calling this function for each, a
fully atom-mapped reaction SMILES can be constructed, allowing for
later reconstruction of atom-mapped geometries.

Heavily based on the implementation in Colin Grambow's `ard_gsm`
package: https://github.com/cgrambow/ard_gsm/tree/v1.0.0
"""
function atom_map_smiles(frame::Dict{String, Any}, smi::String)
    atoms_in_mol_true = Dict{String, Int}()
    for i in 1:frame["N_atoms"]
        elem = frame["arrays"]["species"][i]
        atoms_in_mol_true[elem] = get(atoms_in_mol_true, elem, 0) + 1
    end

    mol_sanitised = rdChem.MolFromSmiles(smi)
    mol_sanitised = rdChem.AddHs(mol_sanitised)
    atoms_in_mol_sani = Dict{String, Int}()
    for atom in mol_sanitised.GetAtoms()
        elem = pyconvert(String, atom.GetSymbol())
        atoms_in_mol_sani[elem] = get(atoms_in_mol_sani, elem, 0) + 1
    end
    if atoms_in_mol_true != atoms_in_mol_sani
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

    return pyconvert(String, rdChem.MolToSmiles(mol_sanitised))
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
function atom_map_frame(am_smi::String, frame::Dict{String, Any})
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
    mol_target_sb = rdChem.Mol(mol_target)
    for bond in mol_target_sb.GetBonds()
        bond.SetBondType(rdChem.rdchem.BondType."SINGLE")
    end
    for atom in mol_target_sb.GetAtoms()
        atom.SetAtomMapNum(0)
    end
    
    match = pyconvert(Vector, mol_target_sb.GetSubstructMatch(mol_template))
    if pyconvert(Int, mol_template.GetNumAtoms()) != length(match)
        println(mol_template.GetNumAtoms())
        throw(ErrorException("Incorrect number of atoms when matching substruct during atom mapping."))
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