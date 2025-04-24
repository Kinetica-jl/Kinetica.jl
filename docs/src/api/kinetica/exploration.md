# Kinetica.jl API

## CRN Representation

### Representing Chemical Species

```@docs
SpeciesData
Base.push!(::SpeciesData, ::String, ::Dict{String, Any}, ::Int)
Base.push!(::SpeciesData, ::Vector{String}, ::Vector{Any}, ::Int)
Base.push!(::SpeciesData, ::String, ::Int)
push_unique!
```

### Representing Reactions

```@docs
RxData
Base.push!(::RxData{iType, fType}, ::SpeciesData, ::Vector{Vector{String}}, ::Vector{Vector{String}}, ::Vector{Dict{String, Any}}, ::Vector{Dict{String, Any}}, ::Vector{fType}, ::Int) where {iType, fType <: AbstractFloat}
Base.splice!(::RxData, ::Vector{Int})
get_rhash
get_reverse_rhash
```

### CRN Initialisation

```@docs
init_network
```

## CDE Interface

```@docs
CDE
ingest_cde_run
import_mechanism
import_mechanism!
import_network
```

## Exploration

```@docs
Kinetica.ExploreLoc
pathof(::Kinetica.ExploreLoc)
Kinetica.find_current_loc
DirectExplore
IterativeExplore
explore_network
Kinetica.load_past_seeds
Kinetica.load_current_seeds
Kinetica.identify_next_seeds
```

## Molecule System

```@docs
system_from_smiles
system_from_mols
Kinetica.molsys_opt
```