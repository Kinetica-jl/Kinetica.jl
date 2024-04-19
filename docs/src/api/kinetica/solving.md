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