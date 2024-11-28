using Kinetica
using Test
using SafeTestsets
using Random
Random.seed!(12345)

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "Main"
    @safetestset "ConditionSet Tests" begin
        include("Main/conditions.jl")
    end
end

if GROUP == "All" || GROUP == "Python"
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
end