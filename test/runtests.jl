using Kinetica
using Test
using SafeTestsets

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
end