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