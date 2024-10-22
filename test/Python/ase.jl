using Test
using Kinetica
using PythonCall

@testset "ASE Conversion Tests" begin
    ethene_xyz = xyz_from_smiles("C=C"; generator=:rdkit, seed=10)
    ethene_frame = xyz_to_frame(ethene_xyz)
    atoms = frame_to_atoms(ethene_frame)
    @test pylen(atoms) == 6

    asecalcemt = pyimport("ase.calculators.emt")
    atoms.calc = asecalcemt.EMT()
    newframe = atoms_to_frame(atoms, pyconvert(Float64, atoms.get_potential_energy()))
    @test Float32(newframe["info"]["energy_ASE"]) ≈ Float32(2.41676)

    @test Float32(Kinetica.imaginary_ve_tol(1e-3)) ≈ Float32(2.0445437750827997)
end

@testset "Basic Species Property Calculation" begin
    xyz = xyz_from_smiles("C[C][CH2-]"; generator=:rdkit, seed=10)
    frame = xyz_to_frame(xyz)
    sd = SpeciesData(["C[C][CH2-]"], [frame])
    sd.cache[:mult] = Dict{Int, Int}()
    sd.cache[:charge] = Dict{Int, Int}()
    sd.cache[:formal_charges] = Dict{Int, Vector{Int}}()
    sd.cache[:initial_magmoms] = Dict{Int, Vector{Float64}}()
    @test Kinetica.get_mult!(sd, 1) == 3
    @test Kinetica.get_charge!(sd, 1) == -1
    @test Kinetica.get_formal_charges!(sd, 1) == [0,0,-1,0,0,0,0,0]
    @test Kinetica.get_initial_magmoms!(sd, 1) == [0.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0]

    reac_magmoms = [2.0, 1.0, 0.0, 0.0, 0.0]
    prod_magmoms = [0.0, 0.0, 0.0, 0.0, 0.0]
    @test_throws "Reactant magmoms cannot be corrected" Kinetica.correct_magmoms_for_mult!(reac_magmoms, prod_magmoms, 1)
    reac_magmoms = [2.0, 0.0, 0.0, 0.0, 0.0]
    prod_magmoms = [1.0, 1.0, 0.0, 0.0, 0.0]
    Kinetica.correct_magmoms_for_mult!(reac_magmoms, prod_magmoms, 1)
    @test reac_magmoms == [0.0, 0.0, 0.0, 0.0, 0.0]
    @test prod_magmoms == [1.0, -1.0, 0.0, 0.0, 0.0]

    @test Kinetica.get_hydrogen_idxs("[C:1](=[C:2]([H:5])[H:6])([H:3])[H:4]") == [[5, 6, 3, 4]]
end

@testset "Calculator Builders" begin
    builders = [EMTBuilder, NWChemDFTBuilder, FHIAimsBuilder]
    btypenames = [
        "<class 'ase.calculators.emt.EMT'>",
        "<class 'ase.calculators.nwchem.NWChem'>",
        "<class 'ase.calculators.aims.Aims'>"
    ]
    mkpath("./species_defaults/defaults_2020/tight")
    for (i, b) in enumerate(builders)
        builder = b()
        calc = builder("./", 1, 0)
        @test pyconvert(String, pystr(pytype(calc))) == btypenames[i]
    end
    rm("./species_defaults", recursive=true)
end

@testset "Geometry Optimisation" begin
    xyz = xyz_from_smiles("CC"; generator=:rdkit, seed=10)
    frame = xyz_to_frame(xyz)
    sd = SpeciesData(["CC"], [frame])
    sd.cache[:mult] = Dict{Int, Int}(1 => 1)
    sd.cache[:charge] = Dict{Int, Int}(1 => 0)
    sd.cache[:formal_charges] = Dict{Int, Vector{Int}}(1 => [0,0,0,0,0,0,0,0])
    sd.cache[:initial_magmoms] = Dict{Int, Vector{Float64}}(1 => [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    builder = EMTBuilder()
    conv = Kinetica.geomopt!(sd, 1, builder)
    @test conv
    @test Float32(sd.xyz[1]["info"]["energy_ASE"]) ≈ Float32(1.46196)

    tframe = deepcopy(frame)
    tframe["arrays"]["pos"] .+= [0.5, 0.75, 0.2]
    Kinetica.kabsch_fit!(frame, tframe)
    @test Float32.(frame["arrays"]["pos"]) ≈ Float32.(tframe["arrays"]["pos"])
end

@testset "ASENEBCalculator" begin
    sd = SpeciesData(["CC", "C=C", "[H][H]"], [
        xyz_to_frame(xyz_from_smiles("CC"; generator=:rdkit, seed=10)),
        xyz_to_frame(xyz_from_smiles("C=C"; generator=:rdkit, seed=10)),
        xyz_to_frame(xyz_from_smiles("[H][H]"; generator=:rdkit, seed=10))
    ])
    h2frame_translated = deepcopy(sd.xyz[3])
    h2frame_translated["arrays"]["pos"] .+= [2.5, 2.5, 2.5]
    prodsys = Kinetica.combine_mols([sd.xyz[2], h2frame_translated])
    rd = RxData(sd, [["CC"], ["C=C", "[H][H]"]], [["C=C", "[H][H]"], ["CC"]],
                [deepcopy(sd.xyz[1]), deepcopy(prodsys)],
                [deepcopy(prodsys), deepcopy(sd.xyz[1])],
                [0.0, 0.0])
    builder = EMTBuilder()
    calc = ASENEBCalculator(builder, "./calc"; ftol=0.1, climb_ftol=0.5, n_images=5,
                            imaginary_freq_tol=0.1)
    setup_network!(sd, rd, calc)
    k = calc(; T=300.0, P=1e5)
    # EMT is too terrible to allow for numerical rate constant
    # comparisons here, so we just have to settle for it calculating
    # without failure.
    @test k isa Vector{Float64}
    @test calculate_entropy_enthalpy(calc, 300.0, 1e5) isa Tuple{Vector{Float64}, Vector{Float64}}
    rm("./calc", recursive=true)
end