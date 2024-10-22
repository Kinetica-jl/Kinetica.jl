# [Kinetica.jl API](@id api_solving)

## Simulation Parameters

```@docs
ODESimulationParams
```

## Kinetic Calculators (Kinetica.jl)

```@docs
allows_continuous
setup_network!
has_conditions
Base.splice!(calc::cType, rids::Vector{Int}) where {cType <: Kinetica.AbstractKineticCalculator}
Base.splice!(rd::RxData, calculator::cType, rids::Vector{Int}) where {cType <: Kinetica.AbstractKineticCalculator}
PrecalculatedArrheniusCalculator
```

The [`ASENEBCalculator`](@ref) ASE-driven NEB-based calculator also implemented in Kinetica.jl is significantly more complex, and has its own API page [here](@ref api_ase_calculator).

## Solvers

```@docs
StaticODESolve
VariableODESolve
solve_network
```

## Reaction Filtering

```@docs
RxFilter
```