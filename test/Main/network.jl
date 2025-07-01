using Test
using Kinetica
using PythonCall
using ExtXYZ

@testset "Species Ingest" begin
    xyz_str = xyz_file_to_str("Main/data/C4H8_rad.xyz")
    @test xyz_str == "12\n\nC     -1.9436231500    0.0789638149    0.4136740755\nC     -0.4607265556   -0.2556714467    0.4006128069\nH     -2.2031727850    0.6762224850    1.2935580259\nH     -2.2262110107    0.6434341526   -0.4806093752\nH     -2.5392750754   -0.8388669281    0.4383778012\nC      0.3972874000    1.0072485740    0.3661253753\nH     -0.2404345661   -0.8798204970   -0.4732563400\nC      1.8802039280    0.6726138230    0.3530673849\nH      0.1540962287    1.5987923724   -0.5241762883\nH      2.1627779651    0.1081400219    1.2473548736\nH      2.1397406454    0.0753498304   -0.5268202045\nH      2.4758439220    1.5904408993    0.3283662969\n"
    smi_list, xyz_list = ingest_xyz_system(xyz_str)
    @test length(smi_list) == 1 && length(xyz_list) == 1
    @test smi_list[1] == "CC=CC"
    @test xyz_list[1]["N_atoms"] == 12
    @test xyz_list[1]["arrays"]["species"] == ["C", "C", "H", "H", "H", "C", "H", "C", "H", "H", "H", "H"]
end

@testset "Reaction Ingest" begin
    rdir = "Main/data/fake_crn/level_001/subspace_001/"
    rsmis, rxyzs, rsys, psmis, pxyzs, psys, dH = ingest_cde_run(rdir, 1)

    @test length(rsmis) == 2 && length(psmis) == 2
    @test rsmis[1][1] == "CCCC"
    @test Set(rsmis[2]) == Set(["C/C=C/C", "[H][H]"])
    @test rsmis[1] == psmis[2] && rsmis[2] == psmis[1]

    @test length(rxyzs) == 2 && length(pxyzs) == 2
    @test rxyzs[1][1]["N_atoms"] == 14
    @test pxyzs[1][1]["N_atoms"] + pxyzs[1][2]["N_atoms"] == 14

    @test length(rsys) == 2 && length(psys) == 2
    @test rsys[1]["N_atoms"] == 14 && psys[1]["N_atoms"] == 14

    @test all(dH .== [-5.0, 5.0])
end

@testset "Network Creation" begin
    sd, rd = init_network()
    @test sd.n == 0
    @test rd.nr == 0

    # Tests both constructors and removal of duplicate species.
    sd = SpeciesData("Main/data/2C_C4H10_H2_C4H8rad.xyz")
    @test sd.n == 4
    @test length(sd.toStr) == 4
    @test length(sd.toInt) == 4
    @test length(sd.xyz) == 4
    @test isnothing(sd.surfdata)
    @test isempty(sd.cache)
    all_species = Set(collect(keys(sd.toInt)))
    @test all_species == Set(["C", "CCCC", "[H][H]", "C/C=C/C"])
    
    rsmis, rxyzs, rsys, psmis, pxyzs, psys, dH = ingest_cde_run("Main/data/fake_crn/level_001/subspace_001/", 1)
    rd = RxData(sd, rsmis, psmis, rsys, psys, dH)
    @test rd.nr == 2
    reac_smis = [[sd.toStr[i] for i in id_reac] for id_reac in rd.id_reacs]
    @test reac_smis[1][1] == "CCCC"
    @test Set(reac_smis[2]) == Set(["C/C=C/C", "[H][H]"])
    @test rd.id_prods[1] == rd.id_reacs[2] && rd.id_prods[2] == rd.id_reacs[1]
    @test rd.stoic_reacs == [[1], [1, 1]]
end

@testset "Single Mechanism Import" begin
    loc = Kinetica.ExploreLoc("Main/data/fake_crn", 1, 1)
    @test pathof(loc) == "Main/data/fake_crn/level_001/subspace_001"
    for _ in 1:3 Kinetica.inc_level!(loc) end
    for _ in 1:5 Kinetica.inc_subspace!(loc) end
    @test pathof(loc) == "Main/data/fake_crn/level_004/subspace_006"

    loc = Kinetica.find_current_loc("Main/data")
    @test loc.level == 0 && loc.subspace == 1
    loc = Kinetica.find_current_loc("Main/data/fake_crn")
    @test loc.level == 1 && loc.subspace == 1

    sd, rd = import_mechanism(loc, 1)
    @test sd.n == 3 && rd.nr == 2

    @test bytes2hex(get_rhash(sd, rd, 1)) == "42196fd8d4700f4a2246e07699e7385cd56a3cb49bb9e9a803f33243dc2ad35c" 
    @test get_reverse_rhash(sd, rd, 1) == rd.rhash[2]
end

@testset "Network Modification" begin
    sd, rd = init_network()
    push!(sd, "Main/data/2C_C4H10_H2_C4H8rad.xyz")
    @test sd.n == 5

    sd, rd = init_network()
    push_unique!(sd, "Main/data/2C_C4H10_H2_C4H8rad.xyz")
    @test sd.n == 4

    rsmis, rxyzs, rsys, psmis, pxyzs, psys, dH = ingest_cde_run("Main/data/fake_crn/level_001/subspace_001/", 1)
    push!(rd, sd, rsmis, psmis, rsys, psys, dH)
    @test rd.nr == 2
    # Test reaction uniqueness detection.
    push!(rd, sd, rsmis, psmis, rsys, psys, dH)
    @test rd.nr == 2
    # Disable uniqueness detection.
    push!(rd, sd, rsmis, psmis, rsys, psys, dH; unique_rxns=false)
    @test rd.nr == 4

    remove_rids = [3, 4]
    splice!(rd, remove_rids)
    @test rd.nr == 2

    Kinetica.populate_sd_cache!(sd)
    @test all([haskey(sd.cache, prop) for prop in [
        :vib_energies,
        :symmetry,
        :mult,
        :charge,
        :formal_charges,
        :geometry,
        :initial_magmoms,
        :ads_xyz
    ]])

    get_species_stats!(sd)
    @test sd.cache[:radii][sd.toInt["C"]] ≈ Float32(2.17376)
    @test sd.cache[:weights][sd.toInt["[H][H]"]] ≈ 2.016
    get_species_stats!(sd; refresh=true)
end

@testset "Network Import" begin
    sd, rd = import_network("Main/data/fake_crn")
    @test sd.n == 4 && rd.nr == 2
    @test haskey(sd.toInt, "N#N")
end