using Kinetica
using Test
using SafeTestsets
using Random
Random.seed!(12345)
Kinetica.np.random.seed(12345)

# Python tests - should come first, as they are most likely to have issues and
# many other bits depend on them down the chain.
@safetestset "OpenBabel Tests" begin
    include("Python/openbabel.jl")
end
@safetestset "RDKit Tests" begin
    include("Python/rdkit.jl")
end
@safetestset "ASE Tests" begin
    include("Python/ase.jl")
end
@safetestset "Surface Tests" begin
    include("Python/surfaces.jl")
end
@safetestset "autodE Tests" begin
    include("Python/autode.jl")
end

# Pure Julia tests.
@safetestset "ConditionSet Tests" begin
    include("Main/conditions.jl")
end