# Kinetica.jl API

## CRN Representation

### Representing Chemical Species

```@docs
SpeciesData
Base.push!(::SpeciesData, ::String, ::Dict{String, Any})
Base.push!(::SpeciesData, ::Vector{String}, ::Vector{Any})
Base.push!(::SpeciesData, ::String)
push_unique!
```

### Representing Reactions

```@docs
RxData
Base.push!(::RxData{iType, fType}, ::SpeciesData, ::Vector{Vector{String}}, ::Vector{Vector{String}}, ::Vector{Dict{String, Any}}, ::Vector{Dict{String, Any}}, ::Vector{fType}) where {iType, fType <: AbstractFloat}
Base.splice!(::RxData, ::Vector{Int})
```

### CRN Initialisation

```@docs
init_network
```

## CDE Interface

```@docs
CDE
Kinetica.ingest_cde_run
```

## Exploration

```@docs
DirectExplore
IterativeExplore
explore_network
```

## Molecule System

```@docs
system_from_smiles
system_from_mols
Kinetica.molsys_opt
```