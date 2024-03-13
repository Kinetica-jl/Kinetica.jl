using Test
using Kinetica
using PythonCall

@testset "RDKit Conversion Tests" begin
    ethene_xyz = xyz_from_smiles("C=C"; generator=:rdkit, seed=10)
    ethene_frame = xyz_to_frame(ethene_xyz)
    @test ethene_frame["N_atoms"] == 6

    f = tempname(); xyz_from_smiles("C=C", f; generator=:rdkit, seed=10)
    ethene_xyz_in = xyz_file_to_str(f)
    ethene_frame_in = xyz_to_frame(ethene_xyz_in)
    @test ethene_frame["N_atoms"] == ethene_frame_in["N_atoms"]
    @test ethene_frame["arrays"]["species"] == ethene_frame_in["arrays"]["species"]
    @test Float32.(ethene_frame["arrays"]["pos"]) ≈ Float32.(ethene_frame_in["arrays"]["pos"])

    rdmol = frame_to_rdkit(ethene_frame; with_coords=true)
    rdChem.SanitizeMol(rdmol)
    @test pyconvert(Int, rdmol.GetNumAtoms()) == 6
    @test pyconvert(Int, rdmol.GetNumBonds()) == 5
    conf = rdmol.GetConformer()
    point = conf.GetAtomPosition(0)
    point_coords = [point.x, point.y, point.z]
    @test pyconvert(Vector{Float32}, point_coords) ≈ Float32.(ethene_frame["arrays"]["pos"][:, 1])
end

@testset "RDKit Atom Mapping Tests" begin
    ethene_xyz = xyz_from_smiles("C=C"; generator=:rdkit, seed=10)
    ethene_frame = xyz_to_frame(ethene_xyz)
    amsmi = atom_map_smiles(ethene_frame, "C=C")
    @test amsmi == "[C:1](=[C:2]([H:5])[H:6])([H:3])[H:4]"

    ethene_frame_mod = xyz_to_frame(ethene_xyz)
    ethene_frame_mod["arrays"]["species"] = reverse(ethene_frame_mod["arrays"]["species"])
    ethene_frame_mod["arrays"]["pos"] = reverse(ethene_frame_mod["arrays"]["pos"]; dims=2)
    amframe = atom_map_frame(amsmi, ethene_frame_mod)
    @test ethene_frame["arrays"]["species"] == amframe["arrays"]["species"]
    # Atom mapping hydrogens does not guarantee that they
    # will have consistent indices, so this doesn't always
    # work.
    # @test Float32.(ethene_frame["arrays"]["pos"]) ≈ -Float32.(amframe["arrays"]["pos"])
end