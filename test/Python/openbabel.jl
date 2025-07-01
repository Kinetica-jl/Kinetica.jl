using Test
using Kinetica
using PythonCall

@testset "OpenBabel Conversion Tests" begin
    ethene_frame = frame_from_smiles("C=C")
    @test ethene_frame["N_atoms"] == 6

    ethene_xyz_out = frame_to_xyz(ethene_frame)
    f1 = tempname(); write(f1, ethene_xyz_out)
    ethene_xyz_in = xyz_file_to_str(f1)
    can1 = String(split(pyconvert(String, pybel.readstring("xyz", ethene_xyz_in).write("can")), '\t')[1])
    @test can1 == "C=C"
    f2 = tempname(); xyz_from_smiles("C=C", f2)
    ethene_xyz_in2 = xyz_file_to_str(f2)
    can2 = String(split(pyconvert(String, pybel.readstring("xyz", ethene_xyz_in2).write("can")), '\t')[1])
    @test can2 == "C=C"

    @test_throws ArgumentError xyz_from_smiles("C=C"; seed=1)
    @test_throws ArgumentError xyz_from_smiles("C=C", f2; seed=1)

    f3 = tempname(); system_from_smiles(["CC", "[H][H]", "[CH2][CH2]"], f3)
    smi_list, xyz_list = ingest_xyz_system(xyz_file_to_str(f3))
    @test issetequal(Set(smi_list), Set(["CC", "[H][H]", "C=C"]))
end