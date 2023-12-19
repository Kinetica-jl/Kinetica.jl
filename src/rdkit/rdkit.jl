function frame_to_rdkit(frame::Dict{String, Any})
    pbmol = pybel.readstring("xyz", frame_to_xyz(frame))
    rdmol = rdChem.rdchem.EditableMol(rdChem.rdchem.Mol())
    for i in 1:frame["N_atoms"]
        rdatom = rdChem.rdchem.Atom(frame["arrays"]["species"][i])
        rdatom.SetAtomMapNum(i)
        rdmol.AddAtom(rdatom)
    end

    for obbond in pybel.ob.OBMolBondIter(pbmol.OBMol)
        a1 = obbond.GetBeginAtom()
        a2 = obbond.GetEndAtom()
        idx1 = a1.GetIdx()
        idx2 = a2.GetIdx()
        rdmol.AddBond(idx1-1, idx2-1, rdChem.rdchem.BondType."SINGLE")
    end

    return rdmol.GetMol()
end


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
        elem = atom.GetSymbol()
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

    match = mol_sani_sb.GetSubstructMatch(mol_with_map)
    if mol_with_map.GetNumAtoms() != length(match)
        println(mol_with_map.GetNumAtoms())
        throw(ErrorException("Incorrect number of atoms when matching substruct during atom mapping."))
    end
    for atom in mol_with_map.GetAtoms()
        idx = match[atom.GetIdx() + 1]
        map_num = atom.GetAtomMapNum()
        mol_sanitised.GetAtomWithIdx(idx).SetAtomMapNum(map_num)
    end

    return rdChem.MolToSmiles(mol_sanitised)
end