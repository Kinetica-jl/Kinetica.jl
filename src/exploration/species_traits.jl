abstract type AbstractSpecies end

struct GasSpecies <: AbstractSpecies end
struct SurfaceSpecies <: AbstractSpecies end

"""
    SpeciesStyle(smi::String)

Trait returning correct subtype of `AbstractSpecies` for a given SMILES.

Useful for automatically dispatching on gas/surface SMILES without
having to sacrifice the simplicity of them being `String`s.

Given a function `myfunc(a::Int, smi::String)`, this can be made to
dispatch on SMILES type by piping it to
`myfunc(SpeciesStyle(smi), a::Int, smi::String)` and creating methods
for `myfunc(GasSpecies, a::Int, smi::String)` and
`myfunc(SurfaceSpecies, a::Int, smi::String)`.
"""
function SpeciesStyle(smi::String)
    if "X_" in smi
        return SurfaceSpecies
    else
        return GasSpecies
    end
end