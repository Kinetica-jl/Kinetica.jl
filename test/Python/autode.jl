using Test
using Kinetica
using PythonCall

@testset "autodE Conversion Tests" begin
    xyz = xyz_from_smiles("C[C][CH2-]"; generator=:rdkit, seed=10)
    frame = xyz_to_frame(xyz)
    ademol = Kinetica.frame_to_autode(frame; mult=3, chg=-1)
    @test pyconvert(Float32, ademol.angle(0, 1, 2).real) ≈ 1.93862

    ademol.mult = 1
    info = Dict("mult" => pyconvert(Int, ademol.mult), "chg" => pyconvert(Int, ademol.charge))
    newframe = Kinetica.autode_to_frame(ademol; info_dict=info)
    @test newframe["info"]["mult"] == 1
end

@testset "autodE Utils" begin
    ethene_xyz = xyz_from_smiles("C=C"; generator=:rdkit, seed=10)
    ethene_frame = xyz_to_frame(ethene_xyz)
    sd = SpeciesData(["C=C"], [ethene_frame])
    sd.cache[:mult] = Dict{Int, Int}(1 => 1)
    sd.cache[:charge] = Dict{Int, Int}(1 => 0)
    graph = Kinetica.autode_get_graph(sd, 1)
    @test pyconvert(Int, graph.order()) == 6

    ethene_xyz2 = xyz_from_smiles("C=C"; generator=:rdkit, seed=11)
    ethene_frame2 = xyz_to_frame(ethene_xyz2)
    push!(sd, "C=C", ethene_frame2)
    sd.cache[:mult][2] = 1; sd.cache[:charge][2] = 0
    graph2 = Kinetica.autode_get_graph(sd, 2)
    @test Kinetica.autode_is_isomorphic(graph, graph2) == true

    ethene_frame2["arrays"]["pos"][:, 5] .+= 10
    push!(sd, "C=C", ethene_frame2)
    sd.cache[:mult][3] = 1; sd.cache[:charge][3] = 0
    graph3 = Kinetica.autode_get_graph(sd, 3)
    @test Kinetica.autode_is_isomorphic(graph, graph3) == false

    @test Kinetica.autode_frame_symmetry(ethene_frame) == (4, 2)
end

@testset "autodE Conformer Generation" begin
    ethene_xyz = xyz_from_smiles("C=C"; generator=:rdkit, seed=10)
    ethene_frame = xyz_to_frame(ethene_xyz)
    sd = SpeciesData(["C=C"], [ethene_frame])
    sd.cache[:mult] = Dict{Int, Int}(1 => 1)
    sd.cache[:charge] = Dict{Int, Int}(1 => 0)
    sd.cache[:symmetry] = Dict{Int, Int}()
    sd.cache[:geometry] = Dict{Int, Int}()
    Kinetica.autode_conformer_search!(sd, 1)
    @test Float32(sd.xyz[1]["info"]["energy"]) ≈ Float32(-170.64979316539365)
    @test sd.cache[:geometry][1] == 2

    ethyne_xyz = xyz_from_smiles("C#C"; generator=:rdkit, seed=10)
    ethyne_frame = xyz_to_frame(ethyne_xyz)
    push!(sd, "C#C", ethyne_frame)
    sd.cache[:mult][2] = 1; sd.cache[:charge][2] = 0
    Kinetica.autode_conformer_search!(sd, 2)
    @test sd.cache[:geometry][2] == 1

    sys = Kinetica.autode_NCI_conformer_search(sd, [1, 2])
    ediff = sys["info"]["energy"] + 312.3 
    @test abs(ediff) < 0.5

    rm("./conformers", recursive=true)
    rm("./.autode_calculations")
end