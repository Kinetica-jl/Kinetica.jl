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
for `myfunc(::GasSpecies, a::Int, smi::String)` and
`myfunc(::SurfaceSpecies, a::Int, smi::String)`.
"""
function SpeciesStyle(smi::String)
    if !isnothing(match(r"(X\d_\d)", smi))
        return SurfaceSpecies()
    else
        return GasSpecies()
    end
end
 


abstract type AbstractXYZ end

struct FreeXYZ <: AbstractXYZ end
struct AdsorbateXYZ <: AbstractXYZ end
struct OnSurfaceXYZ <: AbstractXYZ end

"""
    XYZStyle(frame::Dict{String, Any})

Trait returning correct subtype of `AbstractXYZ` for a given ExtXYZ `frame`.

Detects if a frame has any periodic boundaries enabled, in which case
it is assumed to represent one or more species, on or above a surface.
"""
function XYZStyle(frame::Dict{String, Any})
    if haskey(frame, "pbc") && any(frame["pbc"])
        return OnSurfaceXYZ()
    elseif haskey(frame["info"], "adsorbate") && frame["info"]["adsorbate"] == "true"
        return AdsorbateXYZ()
    else
        return FreeXYZ()
    end
end