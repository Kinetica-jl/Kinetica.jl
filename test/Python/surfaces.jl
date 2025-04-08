using Test
using Kinetica
using PythonCall

@testset "Surface Creation" begin
    surf1 = Surface("Au_fcc111", Kinetica.asebuild.fcc111("Au", (3,3,3)))

    @test surf1.name == "Au_fcc111"
    @test surf1.elements == Set(["Au"])
    sites = keys(atoms_to_frame(surf1.atoms)["info"]["adsorbate_info"]["sites"])
    @test Set(sites) == Set(["ontop", "bridge", "fcc", "hcp"])
    @test Set(values(surf1.sites)) == Set(["ontop", "bridge", "fcc", "hcp"])
    @test all([surf1.sites[v] == k for (k, v) in surf1.siteids])

    atoms = Kinetica.asebuild.fcc111("Ag", (1,1,1))
    frame = atoms_to_frame(atoms)
    ads_info = pop!(frame["info"], "adsorbate_info")
    sites = keys(ads_info["sites"])
    sitedict = Dict{String, Any}()
    for site in sites
        Kinetica.asebuild.add_adsorbate(atoms, "H", 1.0, site)
        sitedict[site] = pyconvert(Vector{Float64}, atoms.get_positions()[-1][0:1])
    end
    surf2 = Surface("Ag_fcc111", frame, sitedict; relative_sites=false)

    @test surf2.name == "Ag_fcc111"
    @test surf2.elements == Set(["Ag"])
    sites = keys(atoms_to_frame(surf1.atoms)["info"]["adsorbate_info"]["sites"])
    @test Set(sites) == Set(["ontop", "bridge", "fcc", "hcp"])
    @test Set(values(surf1.sites)) == Set(["ontop", "bridge", "fcc", "hcp"])
    @test all([surf1.sites[v] == k for (k, v) in surf1.siteids])
end

@testset "SurfaceData Creation" begin
    ase_surfs = [
        Kinetica.asebuild.fcc111("Au", (3,3,3)),
        Kinetica.asebuild.fcc100("Pt", (3,3,3))
    ]
    labels = ["Au_fcc111", "Pt_fcc100"]
    surfaces = [Surface(name, surf) for (name, surf) in zip(labels, ase_surfs)]

    surfdata = SurfaceData(surfaces; sf_samples=100)

    @test surfdata.n == 2
    @test all([surfdata.nameToInt[name] == idx for (idx, name) in enumerate(labels)])
    @test pyconvert(Vector{Set}, surfdata.finder.surface_sites) == [
        Set(["ontop", "fcc", "bridge", "hcp"]),
        Set(["ontop", "hollow", "bridge"])
    ]
end

@testset "ASESurfaceFinder Interface" begin
    ase_surfs = [
        Kinetica.asebuild.fcc100("Au", (3,3,3)),
        Kinetica.asebuild.fcc110("Au", (3,3,3)),
        Kinetica.asebuild.fcc111("Au", (3,3,3)),
    ]
    labels = ["Au_fcc100", "Au_fcc110", "Au_fcc111"]
    surfaces = [Surface(name, surf) for (name, surf) in zip(labels, ase_surfs)]

    surfdata = SurfaceData(surfaces)
    testsys = Kinetica.aseio.read("Python/data/CO+H2O+NHCH3_Au_fcc111_opt.xyz")
    slab, molecules, pred_labels_per_mol = surfdata.finder.predict(testsys)

    @test pylen(slab) == 108
    @test pylen(molecules) == 3
    @test pylen(molecules[0]) == 2
    @test pylen(molecules[1]) == 2
    @test pylen(molecules[2]) == 6
    @test pyconvert(String, pred_labels_per_mol[0][0]["site"]) == "Au_fcc111_ontop"
    @test pyconvert(String, pred_labels_per_mol[1][0]["site"]) == "Au_fcc111_fcc"
    @test pyconvert(String, pred_labels_per_mol[2][0]["site"]) == "Au_fcc111_bridge"
end

@testset "Surface Ingest" begin
    surfdata = SurfaceData([Surface("Au_fcc111", Kinetica.asebuild.fcc111("Au", (3,3,3)))])
    sd, rd = init_network(surfdata)
    loc = Kinetica.ExploreLoc("Python/data/surface_crn", 1, 1)
    import_mechanism!(sd, rd, loc, 1)

    @test sd.n == 2
    @test sd.toStr[1] == "O=C=O"
    @test sd.toStr[2] == "[X1_1][O]=C=O"
    @test sd.xyz[2]["N_atoms"] == 3
    @test rd.nr == 2
    @test rd.mapped_rxns[1] == "[C:1](=[O:2])=[O:3]>>[X1_1]<-[O:2]=[C:1]=[O:3]"
    @test rd.mapped_rxns[2] == "[X1_1]<-[O:2]=[C:1]=[O:3]>>[C:1](=[O:2])=[O:3]"
end

@testset "Surface Adsorption" begin
    surfdata = SurfaceData([Surface("Au_fcc111", Kinetica.asebuild.fcc111("Au", (1,1,3)))])
    sd, rd = init_network(surfdata)
    sd.cache[:symmetry] = Dict{Int, Int}()
    sd.cache[:mult] = Dict{Int, Int}()
    sd.cache[:charge] = Dict{Int, Int}()
    sd.cache[:formal_charges] = Dict{Int, Vector{Int}}()
    sd.cache[:geometry] = Dict{Int, Int}()
    sd.cache[:initial_magmoms] = Dict{Int, Vector{Float64}}()
    sd.cache[:ads_xyz] = Dict{Int, Dict{String, Any}}()

    smis, xyzs = ingest_xyz_system(xyz_file_to_str("Python/data/CO+H2O+NHCH3_Au_fcc111_opt.xyz"), sd.surfdata)
    push_unique!(sd, smis, xyzs)

    # Single molecule, default height.
    frame1 = adsorb_frame(sd, 1)
    @test frame1["N_atoms"] == 29
    @test frame1["arrays"]["pos"][3, 28] - frame1["arrays"]["pos"][3, 3] ≈ 2.12
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
    Kinetica.conformer_search!(sd, 1); Kinetica.conformer_search!(sd, 3)
    Kinetica.get_formal_charges!(sd, 1); Kinetica.get_formal_charges!(sd, 3)
    Kinetica.get_initial_magmoms!(sd, 1); Kinetica.get_initial_magmoms!(sd, 3)
    builder = TBLiteBuilder(; method="GFN1-xTB")
    Kinetica.geomopt!(sd, 1, builder; fmax=1.0); Kinetica.geomopt!(sd, 3, builder; fmax=1.0)
    frame4 = adsorb_two_frames(sd, 1, 3)

    # Two molecules, surface/gas.
    smis, xyzs = ingest_xyz_system(xyz_file_to_str("Python/data/C4H10.xyz"), sd.surfdata)
    push_unique!(sd, smis, xyzs)
    Kinetica.get_mult!(sd, 4); Kinetica.get_charge!(sd, 4)
    Kinetica.conformer_search!(sd, 4); rm("conformers", recursive=true)
    Kinetica.get_formal_charges!(sd, 4); Kinetica.get_initial_magmoms!(sd, 4)
    Kinetica.geomopt!(sd, 4, builder; fmax=1.0)
    frame5 = adsorb_two_frames(sd, 1, 4)
    @test frame5["N_atoms"] == 64
end