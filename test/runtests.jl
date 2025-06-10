using Kinetica
using Test
using SafeTestsets
using Random
Random.seed!(12345)
Kinetica.np.random.seed(12345)

# Network import requires I/O through openbabel and RDKit for 
# SMILES generation, molecule separation and atom mapping, so
# these have to come first.
@safetestset "OpenBabel Tests" begin
    include("Python/openbabel.jl")
end
@safetestset "RDKit Tests" begin
    include("Python/rdkit.jl")
end

@safetestset "Base Reaction Network Tests" begin
    include("Main/network.jl")
end
@safetestset "ConditionSet Tests" begin
    include("Main/conditions.jl")
end

# Further Pythonic functionality.
@safetestset "ASE Tests" begin
    include("Python/ase.jl")
end
@safetestset "Surface Tests" begin
    include("Python/surfaces.jl")
end
@safetestset "autodE Tests" begin
    include("Python/autode.jl")
end

