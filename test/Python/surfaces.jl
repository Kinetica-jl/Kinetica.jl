using Test
using Kinetica
using PythonCall
using ExtXYZ

const SampleBounds = Kinetica.asesf.SampleBounds

@testset "Surface Creation" begin
    surf1 = Surface("Au_fcc111", Kinetica.asebuild.fcc111("Au", (3,3,3));
                    cutoff_mult=1.1, sitecoords=Dict("fcc" => 3, "hcp" => 3),
                    sitebounds=Dict("ontop" => SampleBounds(0.5, 1.5, 2.4)))

    @test surf1.name == "Au_fcc111"
    @test surf1.elements == Set(["Au"])
    sites = keys(atoms_to_frame(surf1.atoms)["info"]["adsorbate_info"]["sites"])
    @test Set(sites) == Set(["ontop", "bridge", "fcc", "hcp"])
    @test Set(values(surf1.sites)) == Set(["ontop", "bridge", "fcc", "hcp"])
    @test all([surf1.sites[v] == k for (k, v) in surf1.siteids])
    @test surf1.sitecoords == Dict("ontop" => -1, "bridge" => -1, "fcc" => 3, "hcp" => 3)
    @test pyconvert(Bool, surf1.sitebounds["ontop"].is_ovoid) == false
    @test all(isnothing.([surf1.sitebounds["fcc"], surf1.sitebounds["hcp"], surf1.sitebounds["bridge"]]))

    atoms = Kinetica.asebuild.fcc111("Ag", (1,1,1))
    frame = atoms_to_frame(atoms)
    ads_info = pop!(frame["info"], "adsorbate_info")
    sites = keys(ads_info["sites"])
    sitedict = Dict{String, Any}()
    for site in sites
        Kinetica.asebuild.add_adsorbate(atoms, "H", 1.0, site)
        sitedict[site] = pyconvert(Vector{Float64}, atoms.get_positions()[-1][0:1])
    end
    surf2 = Surface("Ag_fcc111", frame, sitedict; 
                    relative_sites=false, sitecoords=Dict("bridge" => 2),
                    sitebounds=Dict("bridge" => SampleBounds(0.3, z_min=-0.1, z_mid=0.8, z_max=1.4)))

    @test surf2.name == "Ag_fcc111"
    @test surf2.elements == Set(["Ag"])
    sites = keys(atoms_to_frame(surf2.atoms)["info"]["adsorbate_info"]["sites"])
    @test Set(sites) == Set(["ontop", "bridge", "fcc", "hcp"])
    @test Set(values(surf2.sites)) == Set(["ontop", "bridge", "fcc", "hcp"])
    @test all([surf2.sites[v] == k for (k, v) in surf2.siteids])
    @test surf2.sitecoords == Dict("ontop" => -1, "bridge" => 2, "fcc" => -1, "hcp" => -1)
    @test pyconvert(Bool, surf2.sitebounds["bridge"].is_ovoid)
    @test all(isnothing.([surf2.sitebounds["fcc"], surf2.sitebounds["hcp"], surf2.sitebounds["ontop"]]))
end

@testset "SurfaceData Creation" begin
    surfaces = [
        Surface("Au_fcc111", 
            read_frame("Python/data/Au_fcc111_113_emt.xyz"), 
            Dict("ontop" => [0.0, 0.0],"bridge" => [0.5, 0.0],"fcc" => [1/3, 1/3],"hcp" => [2/3, 2/3]),
            sitecoords=Dict("fcc" => 3, "hcp" => 3),
            sitebounds=Dict("ontop" => SampleBounds(0.5, 1.5, 2.4))
        ),
        Surface("Pt_fcc100", 
            Kinetica.asebuild.fcc100("Pt", (3,3,3)), 
            cutoff_mult=1.2,
            sitecoords=Dict("hollow" => 4),
            sitebounds=Dict(
                "hollow" => SampleBounds(0.5, z_min=0.4, z_mid=0.8, z_max=1.75),
                "bridge" => SampleBounds(0.45, z_min=0.75, z_mid=1.4, z_max=2.0)
            )
        )
    ]
    surfdata = SurfaceData(surfaces; sf_samples=100, sample_defaults=SampleBounds(0.4, 1.5, 2.4))

    @test surfdata.n == 2
    @test all([surfdata.nameToInt[name] == idx for (idx, name) in enumerate(["Au_fcc111", "Pt_fcc100"])])
    @test pyconvert(Vector{Set}, surfdata.finder.surface_sites) == [
        Set(["ontop", "fcc", "bridge", "hcp"]),
        Set(["ontop", "hollow", "bridge"])
    ]
    @test surfdata.surfaces[1].sitecoords == Dict(
        "ontop" => 1,
        "fcc" => 3,
        "hcp" => 3,
        "bridge" => 2
    )
    @test pyconvert(Bool, surfdata.surfaces[1].sitebounds["bridge"].r_max == 0.4)
    @test pyconvert(Bool, surfdata.surfaces[1].sitebounds["ontop"].r_max == 0.5)
    @test surfdata.surfaces[2].sitecoords == Dict(
        "ontop" => 1,
        "bridge" => 2,
        "hollow" => 4
    )
    @test pyconvert(Bool, surfdata.surfaces[2].sitebounds["ontop"].r_max == 0.4)
    @test pyconvert(Bool, surfdata.surfaces[2].sitebounds["hollow"].is_ovoid)
    @test pyconvert(Bool, surfdata.surfaces[2].sitebounds["bridge"].is_ovoid)
end

@testset "ASESurfaceFinder Interface" begin
    ase_surfs = [
        Kinetica.asebuild.fcc100("Au", (1,1,3)),
        Kinetica.asebuild.fcc110("Au", (1,1,3)),
        Kinetica.asebuild.fcc111("Au", (1,1,3)),
    ]
    labels = ["Au_fcc100", "Au_fcc110", "Au_fcc111"]
    surfaces = [Surface(name, surf; cutoff_mult=1.1) for (name, surf) in zip(labels, ase_surfs)]

    surfdata = SurfaceData(surfaces)
    testsys = Kinetica.aseio.read("Python/data/CO+H2O+NHCH3_Au_fcc111_opt.xyz")
    slab, molecules, pred_labels_per_mol = Kinetica.predict_surface_sites(surfdata, testsys)

    @test pylen(slab) == 108
    @test pylen(molecules) == 3
    @test pylen(molecules[0]) == 2
    @test pylen(molecules[1]) == 2
    @test pylen(molecules[2]) == 6
    @test pyconvert(String, pred_labels_per_mol[0][0]["site"]) == "Au_fcc111_ontop"
    @test pyconvert(String, pred_labels_per_mol[1][0]["site"]) == "Au_fcc111_fcc"
    @test pyconvert(String, pred_labels_per_mol[2][0]["site"]) == "Au_fcc111_bridge"
    @test pyconvert(Int, pred_labels_per_mol[0][0]["coordination"]) == 1
    @test pyconvert(Int, pred_labels_per_mol[1][0]["coordination"]) == 3
    @test pyconvert(Int, pred_labels_per_mol[2][0]["coordination"]) == 2
    @test pyconvert(Float32, pred_labels_per_mol[0][0]["height"]) ≈ 2.42497
    @test pyconvert(Float32, pred_labels_per_mol[1][0]["height"]) ≈ 1.77101
    @test pyconvert(Float32, pred_labels_per_mol[2][0]["height"]) ≈ 1.91647
end

@testset "Surface Ingest" begin
    surfdata = SurfaceData([
        Surface("Au_fcc111", Kinetica.asebuild.fcc111("Au", (3,3,3));
                sitebounds=Dict(
                    "ontop"=>SampleBounds(0.5, 1.5, 2.4)
                ))
    ])
    sd, rd = init_network(surfdata)
    loc = Kinetica.ExploreLoc("Python/data/surface_crn", 1, 1)
    import_mechanism!(sd, rd, loc, 1)

    @test sd.n == 4
    @test sd.toStr[1] == "O=C=O"
    @test sd.toStr[2] == "[X1_4][H]"
    @test sd.toStr[3] == "[X1_4][O]=C=O"
    @test sd.toStr[4] == "[H]"
    @test sd.xyz[3]["N_atoms"] == 3
    @test rd.nr == 4
    @test rd.mapped_rxns[1] == "[O:2]=[C:1]=[O:3]>>[X1_4][O:3]=[C:1]=[O:2]"
    @test rd.mapped_rxns[2] == "[X1_4][H:1]>>[H:1]"
end

@testset "Surface Adsorption" begin
    surf_atoms = Kinetica.asebuild.fcc111("Au", (1,1,3))
    c = Kinetica.aseconstraints.FixAtoms(indices=[0, 1])
    surf_atoms.set_constraint(c)
    surf = Surface("Au_fcc111", surf_atoms, sitebounds=Dict(
        "ontop" => SampleBounds(0.6, 1.5, 2.4),
        "bridge" => SampleBounds(0.35, z_min=0.7, z_mid=1.6, z_max=2.2),
        "fcc" => SampleBounds(0.35, z_min=0.4, z_mid=1.4, z_max=2.0),
        "hcp" => SampleBounds(0.35, z_min=0.4, z_mid=1.4, z_max=2.0)
    ))
    surfdata = SurfaceData([surf], sf_samples=1000)
    sd, rd = init_network(surfdata)
    Kinetica.populate_sd_cache!(sd)

    smis, xyzs = ingest_xyz_system(xyz_file_to_str("Python/data/CO+H2O+NHCH3_Au_fcc111_opt.xyz"), sd.surfdata)
    push_unique!(sd, smis, xyzs)

    # Single molecule, default height.
    frame1 = adsorb_frame(sd, 1)
    @test frame1["N_atoms"] == 29
    @test frame1["arrays"]["pos"][3, 28] - frame1["arrays"]["pos"][3, 3] ≈ 2.12
    @test frame1["arrays"]["fixed_pos"] == vcat(repeat([1, 1, 0], 9), [0, 0])
    # Single molecule, custom height.
    frame2 = adsorb_frame(sd, 1, [1.8])
    @test frame2["arrays"]["pos"][3, 28] - frame2["arrays"]["pos"][3, 3] ≈ 1.8
    # Single molecule, exceeds bounds of unit cell.
    frame3 = adsorb_frame(sd, 3)
    @test frame3["N_atoms"] == 33

    # Two molecules, surface/surface.
    # Adsorbates need to be optimised or placement on surface will be wrong.
    Kinetica.get_mult!(sd, 1); Kinetica.get_charge!(sd, 1)
    Kinetica.get_mult!(sd, 3); Kinetica.get_charge!(sd, 3)
    Kinetica.conformer_search!(sd, 1; n_samples=13); Kinetica.conformer_search!(sd, 3; n_samples=13)
    Kinetica.get_formal_charges!(sd, 1); Kinetica.get_formal_charges!(sd, 3)
    Kinetica.get_initial_magmoms!(sd, 1); Kinetica.get_initial_magmoms!(sd, 3)
    builder = TBLiteBuilder(; method="GFN1-xTB")
    Kinetica.geomopt!(sd, 1, builder; fmax=1.0); Kinetica.geomopt!(sd, 3, builder; fmax=1.0)
    frame4 = adsorb_two_frames(sd, 1, 3)
    @test frame4["N_atoms"] == 56
    @test frame4["info"]["unit_cell_mult"] == 4
    @test frame4["arrays"]["pos"][:, 49] ≈ Float32[0.0, 0.0, 16.83118]
    @test frame4["arrays"]["pos"][:, 51] ≈ Float32[5.76999, 2.49847, 16.781178]
    @test frame4["arrays"]["fixed_pos"] == vcat(repeat([1, 1, 0], 16), zeros(Int, 8))
    @test frame4["info"]["ads_sitetags"] == ["X1_4->49", "X1_1->51"]
    @test frame4["info"]["n_adsorbates"] == 2
    @test frame4["info"]["n_gas_species"] == 0

    # Two molecules, surface/gas.
    smis, xyzs = ingest_xyz_system(xyz_file_to_str("Python/data/C4H10.xyz"), sd.surfdata)
    push_unique!(sd, smis, xyzs)
    Kinetica.get_mult!(sd, 4); Kinetica.get_charge!(sd, 4)
    Kinetica.conformer_search!(sd, 4); rm("conformers", recursive=true)
    Kinetica.get_formal_charges!(sd, 4); Kinetica.get_initial_magmoms!(sd, 4)
    Kinetica.geomopt!(sd, 4, builder; fmax=1.0)
    frame5 = adsorb_two_frames(sd, 1, 4)
    @test frame5["N_atoms"] == 64
    @test frame5["arrays"]["pos"][:, 49] ≈ Float32[0.0, 0.0, 16.83118]
    @test frame5["info"]["ads_sitetags"] == ["X1_4->49"]
    @test frame5["info"]["n_adsorbates"] == 1
    @test frame5["info"]["n_gas_species"] == 1
    # Test lowest atom height of gas adsorbate as autodE conformer search 
    # can't be trusted to place it deterministically.
    gas_idxs = collect(51:64)
    surf_idxs = findall(x -> x == "Au", frame5["arrays"]["species"])
    @test minimum(frame5["arrays"]["pos"][3, gas_idxs]) - maximum(frame5["arrays"]["pos"][3, surf_idxs]) >= 5.0
    @test frame5["arrays"]["fixed_pos"] == vcat(repeat([1, 1, 0], 16), zeros(Int, 16))
end

@testset "Surface Reaction Endpoint Matching" begin
    surf_atoms = Kinetica.asebuild.fcc100("Pt", (1,1,3))
    c = Kinetica.aseconstraints.FixAtoms(indices=[0, 1])
    surf_atoms.set_constraint(c)
    surfdata = SurfaceData([Surface("Pt_fcc100", surf_atoms, cutoff_mult=1.2)])
    sd, rd = init_network(surfdata)
    Kinetica.populate_sd_cache!(sd)
    smis, xyzs = ingest_xyz_system(xyz_file_to_str("Python/data/CH3CH2NH2.xyz"), sd.surfdata)
    push_unique!(sd, smis, xyzs)
    smis, xyzs = ingest_xyz_system(xyz_file_to_str("Python/data/CH3CH2NH+H_PtFCC100.xyz"), sd.surfdata)
    push_unique!(sd, smis, xyzs)
    smis, xyzs = ingest_xyz_system(xyz_file_to_str("Python/data/Hgas_PtFCC100.xyz"), sd.surfdata)
    push_unique!(sd, smis, xyzs)

    builder = TBLiteBuilder(; method="GFN1-xTB")
    conv = [false for _ in 1:sd.n]
    for i in 1:sd.n
        Kinetica.get_mult!(sd, i)
        Kinetica.get_charge!(sd, i)
        # Lift degeneracy of symmetric surface rotamers to be deterministic.
        Kinetica.conformer_search!(sd, i; n_samples=13)
        Kinetica.get_formal_charges!(sd, i)
        Kinetica.get_initial_magmoms!(sd, i)
        conv[i] = Kinetica.geomopt!(sd, i, builder; fmax=1.0)
    end
    @test all(conv)

    # Expansion of ads/gas surface to match larger ads/ads surface
    prod1 = adsorb_two_frames(sd, 2, 3)
    reac1 = adsorb_two_frames(sd, 2, 4)
    @test reac1["N_atoms"] > prod1["N_atoms"]
    Kinetica.scale_surface_to_match!(prod1, reac1, sd.surfdata.surfaces[get_surfid(sd.toStr[2])])
    @test reac1["N_atoms"] == prod1["N_atoms"]

    # Placement of expanded surface under gas reactant to match ads/ads surface
    reac2 = deepcopy(sd.xyz[1])
    Kinetica.add_surface_beneath!(reac2, sd.surfdata.surfaces[get_surfid(sd.toStr[2])], 3)
    @test_throws "already contains a surface" Kinetica.add_surface_beneath!(reac2, sd.surfdata.surfaces[get_surfid(sd.toStr[2])], 3)
    @test reac2["N_atoms"] < prod1["N_atoms"]
    Kinetica.scale_surface_to_match!(reac2, prod1, sd.surfdata.surfaces[get_surfid(sd.toStr[2])])
    @test reac2["N_atoms"] == prod1["N_atoms"]

    # Correct surface matching from creation of gas reactant on surface
    reac3 = deepcopy(sd.xyz[1])
    Kinetica.add_surface_beneath!(reac3, sd.surfdata.surfaces[get_surfid(sd.toStr[2])], prod1["info"]["unit_cell_mult"])
    @test reac3["N_atoms"] == prod1["N_atoms"]
end

@testset "Surface Atom Mapping" begin
    surfdata = SurfaceData([Surface("Pt_fcc100", Kinetica.asebuild.fcc100("Pt", (1,1,3)))])
    sd, rd = init_network(surfdata)
    Kinetica.populate_sd_cache!(sd)
    smis, xyzs = ingest_xyz_system(xyz_file_to_str("Python/data/CH3CH2NH+H_PtFCC100.xyz"), sd.surfdata)
    push_unique!(sd, smis, xyzs)
    smis, xyzs = ingest_xyz_system(xyz_file_to_str("Python/data/CH3CH2NH2.xyz"), sd.surfdata)
    push_unique!(sd, smis, xyzs)

    builder = TBLiteBuilder(; method="GFN1-xTB")
    for i in 1:sd.n
        Kinetica.get_mult!(sd, i)
        Kinetica.get_charge!(sd, i)
        Kinetica.conformer_search!(sd, i; n_samples=13)
        Kinetica.get_formal_charges!(sd, i)
        Kinetica.get_initial_magmoms!(sd, i)
        Kinetica.geomopt!(sd, i, builder; fmax=1.0)
    end

    amsmi = atom_map_smiles(sd.xyz[1], sd.toStr[1])
    @test amsmi == "[C:1]([C:2]([N:3]([X1_3])[H:9])([H:7])[H:8])([H:6])([H:5])[H:4]"
    smi_no_surf = "CCN"
    @test_throws "gas-phase SMILES from a surface-bound geometry" atom_map_smiles(sd.xyz[1], smi_no_surf)
    @test_throws "gas-phase SMILES from a surface-bound geometry" atom_map_smiles(sd.cache[:ads_xyz][1], smi_no_surf)
    @test_throws "mapped surface SMILES from a surface-bound species" atom_map_smiles(sd.cache[:ads_xyz][1], sd.toStr[1])

    amsmi2 = "[X1_3][N:1]([C:3]([C:2]([H:9])([H:8])[H:7])([H:6])[H:5])[H:4]"
    amsmi2_no_surf = "[N:1]([C:3]([C:2]([H:9])([H:8])[H:7])([H:6])[H:5])[H:4]"
    amframe = atom_map_frame(amsmi2, sd.xyz[1])
    @test amframe["N_atoms"] == 9
    @test amframe["arrays"]["species"] == ["N", "C", "C", "H", "H", "H", "H", "H", "H"]
    @test amframe["info"]["adsorbate"] == "true"
    amframe_with_surf = atom_map_frame(amsmi2, sd.cache[:ads_xyz][1])
    @test amframe_with_surf["N_atoms"] == 57
    @test amframe_with_surf["arrays"]["species"] == vcat(["N", "C", "C", "H", "H", "H", "H", "H", "H"], ["Pt" for _ in 1:48])
    @test_throws "adsorbate geometry from a gas-phase SMILES" atom_map_frame(amsmi2_no_surf, sd.xyz[1])
    @test_throws "gas-phase geometry from a surface-bound SMILES" atom_map_frame(amsmi2, sd.xyz[3])

    # Multiply-adsorbed surfaces also need testing since they sometimes
    # exhibit strange behaviour.
    adssys = adsorb_two_frames(sd, 1, 2)
    adssys_smi = join([sd.toStr[1], sd.toStr[2]], ".")
    @test_throws "Remove surface atoms from geometry" atom_map_smiles(adssys, adssys_smi)
    adssys_nosurf = deepcopy(adssys)
    Kinetica.remove_surface_atoms!(adssys_nosurf, sd.surfdata, get_surfid(adssys_smi), true)
    adssys_amsmi = atom_map_smiles(adssys_nosurf, adssys_smi)
    @test adssys_amsmi == "[C:1]([C:2]([N:3]([X1_3])[H:9])([H:7])[H:8])([H:6])([H:5])[H:4].[X1_3][H:10]"

    adssys_amsmi2 = "[X1_3][H:1].[X1_3][N:4]([C:3]([C:2]([H:5])([H:6])[H:7])([H:8])[H:9])[H:10]"
    adssys2 = atom_map_frame(adssys_amsmi2, adssys)
    @test adssys["arrays"]["pos"][:, 37] == adssys2["arrays"]["pos"][:, 1] # H:1 moved

    adssys_amsmi3 = "[X1_1][H:1].[X1_1][N:4]([C:3]([C:2]([H:7])([H:6])[H:5])([H:8])[H:9])[H:10]" # H:5 and H:7 switched
    adssys3 = atom_map_frame(adssys_amsmi3, adssys)
    hidxs = Kinetica.get_hydrogen_idxs(adssys_amsmi3)
    Kinetica.permute_hydrogens!(adssys3, hidxs, adssys2)
    @test all(adssys3["arrays"]["pos"][:, 5] .≈ adssys2["arrays"]["pos"][:, 5]) # should be equivalent after swap.
end 

@testset "Full Surface CRN Calculation" begin
    surf_frame = read_frame("Python/data/Au_fcc111_113_emt.xyz")
    surfsites = Dict(
        "ontop" => [0.0, 0.0],
        "bridge" => [0.5, 0.0],
        "fcc" => [1/3, 1/3],
        "hcp" => [2/3, 2/3]
    )
    bounds = Dict(
        "ontop" => SampleBounds(0.5, 1.5, 2.4),
        "bridge" => SampleBounds(0.35, z_min=0.7, z_mid=1.6, z_max=2.2),
        "fcc" => SampleBounds(0.35, z_min=0.4, z_mid=1.4, z_max=2.0),
        "hcp" => SampleBounds(0.35, z_min=0.4, z_mid=1.4, z_max=2.0)
    )
    coords = Dict(
        "fcc" => 3,
        "hcp" => 3,
        "bridge" => 2
    )
    surf = Surface("Au_fcc111", surf_frame, surfsites; 
                   sitebounds=bounds, sitecoords=coords, cutoff_mult=1.2)
    surfdata = SurfaceData([surf])
    sd, rd = init_network(surfdata)
    loc = Kinetica.ExploreLoc("Python/data/surface_crn", 1, 1)
    import_mechanism!(sd, rd, loc, 1)

    builder = EMTBuilder()
    calcdir_head = "./calc_surface_tests"
    calc = ASENEBCalculator(builder, calcdir_head; n_images=5, ftol=0.01, climb=true, neb_optimiser="fire")
    setup_network!(sd, rd, calc) 

    rates = calc(; T=300.0, P=1e5)
    @test length(rates) == 4
    # EMT not reliable enough for rate tests, as long as we can get some
    # values, it's good enough for now.
    # @test all(abs.(log10.(rates) .- [12.0, -40.0, 20.0, 7.0]) .< 1.0)

    rm("./calc_surface_tests", recursive=true)
end