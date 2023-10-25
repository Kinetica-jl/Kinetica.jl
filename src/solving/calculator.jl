abstract type AbstractKineticCalculator end

"""
    allows_continuous(calculator)

Indicates whether a `calculator` is allowed to function in continuous variable rate constant simulations.

Should only ever return `true` for calculators with analytic
expressions for rate constants, as only these can be encoded
within Symbolics and MTK.

Calculators relying on external function calls should therefore
return `false`, as these calls cannot be encoded. These calculators
will only be functional for discrete rate constant update simulations.
"""
function allows_continuous end


"""
    setup_network!(sd, rd, calculator)

Sets up a network for calculation with the provided `calculator`.

Must be implemented for each `calculator`. Called at the
start of a kinetic simulation to check compatibility of
network with calculator, and to populate network-dependent
fields within the calculator.
"""
function setup_network! end


"""
    splice!(calculator, rids)

Removes all information from the reactions at `rids` from `calculator`.

Useful in conjunction with `splice!(rd, rids)` for removing
e.g. low-rate reactions that have been removed from a network
from that network's calculator.

Relies on a calculator-specific implementation, as it
directly modifies fields of that calculator.
"""
function Base.splice!(calc::cType, rids::Vector{Int}) where {cType <: AbstractKineticCalculator}
    splice!(calc, rids)
end


"""
    splice!(rd, calculator, rids)

Convenience wrapper for deleting reaction data from both a network and its calculator.
"""
function Base.splice!(rd::RxData, calculator::cType, rids::Vector{Int}) where {cType <: AbstractKineticCalculator}
    splice!(rd, rids)
    splice!(calculator, rids)
end


"""
Placeholder kinetic calculator.

Kinetic calculator that always returns an array of original
rates when called. Not intended for use in calculations, mainly
for debugging and reference.

Implemented conditions:
* Temperature (`T`, unit: K)
* Volume (`V`, unit: dm^3)

Has support for dispatching with/without a maximum rate constant
`k_max` and scaling by time unit `t_unit` (assuming rates are
provided in units of /s).
"""
mutable struct DummyKineticCalculator{kmType, uType, tType} <: AbstractKineticCalculator
    rates::Vector{uType}
    k_max::kmType
    t_unit::String
    t_mult::tType
end

"""
    calculator = DummyKineticCalculator(rates[, k_max, t_unit])

Outer constructor method for placeholder kinetic calculator.
"""
function DummyKineticCalculator(rates::Vector{uType}; 
        k_max::Union{Nothing, uType}=nothing, t_unit::String="s") where {uType <: AbstractFloat}

    t_mult = tconvert(t_unit, "s")
    return DummyKineticCalculator(rates, k_max, t_unit, t_mult)
end

function setup_network!(sd::SpeciesData, rd::RxData, calc::DummyKineticCalculator)
    if length(calc.rates) != rd.nr
        throw(ArgumentError("Number of rates ($(length(calc.rates))) does not match number of reactions in `RxData` ($(rd.nr))"))
    end
end

function Base.splice!(calc::DummyKineticCalculator, rids)
    splice!(calc.rates, rids)
end

# Dispatched with k_max awareness.
"""
    rates = calculator(; T, V)

Calculate rates with placeholder kinetic calculator.

Accepts temperature (`T`) and volume (`V`) as keyword arguments.
Has methods implemeted for either or both being present, but will
fail to calculate if neither are present.

Automatically dispatches to a method with correct formula for
`k_max`-aware calculation if this is defined in the underlying
`DummyKineticCalculator`.
"""
function (calc::DummyKineticCalculator{uType, uType, tType})(; T=nothing, V=nothing) where {uType, tType}
    return calc(T, V)
end
function (calc::DummyKineticCalculator{uType, uType, tType})(T::Number, V::Number) where {uType, tType}
    return 1.0 ./ ((1.0 / calc.k_max) .+ (1.0 ./ calc.rates)) * calc.t_mult
end
function (calc::DummyKineticCalculator{uType, uType, tType})(T::Number, ::Nothing) where {uType, tType}
    return calc(T, 0.0)
end
function (calc::DummyKineticCalculator{uType, uType, tType})(::Nothing, V::Number) where {uType, tType}
    return calc(0.0, V)
end

# Dispatched without k_max awareness.
function (calc::DummyKineticCalculator{Nothing, uType, tType})(; T=nothing, V=nothing) where {uType, tType}
    return calc(T, V)
end
function (calc::DummyKineticCalculator{Nothing, uType, tType})(T::Number, V::Number) where {uType, tType}
    return calc.rates * calc.t_mult
end
function (calc::DummyKineticCalculator{Nothing, uType, tType})(T::Number, ::Nothing) where {uType, tType}
    return calc(T, 0.0)
end
function (calc::DummyKineticCalculator{Nothing, uType, tType})(::Nothing, V::Number) where {uType, tType}
    return calc(0.0, V)
end

function has_conditions(::DummyKineticCalculator, symbols::Vector{Symbol})
    return all([sym in [:T, :V] for sym in symbols])
end

allows_continuous(::DummyKineticCalculator) = true


"""
Arrhenius theory kinetic calculator for precalculated reactions.

Kinetic calculator that uses the Arrhenius equation to
determine rates of reaction. Requires prior specification of
reaction activation energies (`Ea`) and Arrhenius prefactors
(`A`), will not calculate/predict these internally.

Implemented conditions:
* Temperature (`T`, unit: K)

Requires:
* Activation energies (`Ea`, unit: J/mol)
* Arrhenius prefactors (`A`, unit: mol dm^-3 s^-1 assuming bimolecular reactions)

Has support for dispatching with/without a maximum rate constant
`k_max` and scaling by time unit `t_unit` (assuming rates are
provided in units of /s).
"""
mutable struct PrecalculatedArrheniusCalculator{kmType, uType, tType} <: AbstractKineticCalculator
    Ea::Vector{uType}
    A::Vector{uType}
    k_max::kmType
    t_unit::String
    t_mult::tType
end

"""
    calculator = PrecalculatedArrheniusCalculator(Ea, A[, k_max, t_unit])

Outer constructor method for Arrhenius theory kinetic calculator.
"""
function PrecalculatedArrheniusCalculator(Ea::Vector{uType}, A::Vector{uType}; 
        k_max::Union{Nothing, uType}=nothing, t_unit::String="s") where {uType <: AbstractFloat}
    
    t_mult = tconvert(t_unit, "s")
    return PrecalculatedArrheniusCalculator(Ea, A, k_max, t_unit, t_mult)
end

function setup_network!(sd::SpeciesData, rd::RxData, calc::PrecalculatedArrheniusCalculator)
    if length(calc.Ea) != rd.nr || length(calc.A) != rd.nr
        throw(ArgumentError("Number of parameters (Ea: $(length(calc.Ea)), A: $(length(calc.A))) does not match number of reactions in `RxData` ($(rd.nr))"))
    end
end

function Base.splice!(calc::PrecalculatedArrheniusCalculator, rids::Vector{Int})
    splice!(calc.Ea, rids)
    splice!(calc.A, rids)
end

# Dispatched with k_max awareness.
"""
    rates = calculator(; T)

Calculate rates with precalculated Arrhenius theory kinetic calculator.

Requires temperature (`T`) as a keyword argument.

Automatically dispatches to a method with correct formula for
`k_max`-aware calculation if this is defined in the underlying
`PrecalculatedArrheniusCalculator`.
"""
function (calc::PrecalculatedArrheniusCalculator{uType, uType, tType})(; T::Number) where {uType, tType}
    k_r = calc.A .* exp.(-calc.Ea / (Constants.R * T)) * Constants.N_A * calc.t_mult
    return 1.0 ./ ((1.0 / calc.k_max) .+ (1.0 ./ k_r))
end

# Dispatched without k_max awareness.
function (calc::PrecalculatedArrheniusCalculator{Nothing, uType, tType})(; T::Number) where {uType, tType}
    k_r = calc.A .* exp.(-calc.Ea / (Constants.R * T)) * Constants.N_A * calc.t_mult
    return k_r
end

function has_conditions(::PrecalculatedArrheniusCalculator, symbols::Vector{Symbol})
    return all([sym in [:T] for sym in symbols])
end

allows_continuous(::PrecalculatedArrheniusCalculator) = true


"""
Pressure-dependent Arrhenius theory kinetic calculator for precalculated reactions.

Kinetic calculator that uses the Arrhenius equation modified
with Lindemann's theory to determine pressure-dependent rates
of reaction. Requires prior specification of reaction 
activation energies (`Ea`) and Arrhenius prefactors (`A`), 
will not calculate/predict these internally.

Implemented conditions:
* Temperature (`T`, unit: K)
* Pressure (`P`, unit: Pa)

Requires:
* Activation energies (`Ea`, unit: J/mol)
* Arrhenius prefactors (`A`, unit: mol dm^-3 s^-1 assuming bimolecular reactions)

Has support for dispatching with/without a maximum rate constant
`k_max` and scaling by time unit `t_unit` (assuming rates are
provided in units of /s).
"""
mutable struct PrecalculatedLindemannCalculator{kmType, uType, tType} <: AbstractKineticCalculator
    Ea::Vector{uType}
    A_0::Vector{uType}
    A_inf::Vector{uType}
    k_max::kmType
    t_unit::String
    t_mult::tType
end

"""
    calculator = PrecalculatedLindemannCalculator(Ea, A[, k_max, t_unit])

Outer constructor method for pressure-dependent Arrhenius theory kinetic calculator.
"""
function PrecalculatedLindemannCalculator(Ea::Vector{uType}, A_0::Vector{uType}, A_inf::Vector{uType}; 
        k_max::Union{Nothing, uType}=nothing, t_unit::String="s") where {uType <: AbstractFloat}
    
    t_mult = tconvert(t_unit, "s")
    return PrecalculatedLindemannCalculator(Ea, A_0, A_inf, k_max, t_unit, t_mult)
end

function setup_network!(sd::SpeciesData, rd::RxData, calc::PrecalculatedLindemannCalculator)
    if any([length(p) != rd.nr for p in [calc.Ea, calc.A_0, calc.A_inf]])
        throw(ArgumentError("Number of parameters (Ea: $(length(calc.Ea)), A_0: $(length(calc.A_0)), A_inf: $(length(calc.A_inf))) does not match number of reactions in `RxData` ($(rd.nr))"))
    end
end

function Base.splice!(calc::PrecalculatedLindemannCalculator, rids::Vector{Int})
    splice!(calc.Ea, rids)
    splice!(calc.A_0, rids)
    splice!(calc.A_inf, rids)
end

# Dispatched with k_max awareness.
"""
    rates = calculator(; T, P)

Calculate rates with pressure-dependent Arrhenius theory kinetic calculator.

Requires temperature (`T`) and pressure (`P`) as keyword arguments.

Automatically dispatches to a method with correct formula for
`k_max`-aware calculation if this is defined in the underlying
`PrecalculatedLindemannCalculator`.
"""
function (calc::PrecalculatedLindemannCalculator{uType, uType, tType})(; T::Number, P::Number) where {uType, tType}
    throw(ErrorException("Lindemann rate constants not implemented yet."))
end

# Dispatched without k_max awareness.
function (calc::PrecalculatedLindemannCalculator{Nothing, uType, tType})(; T::Number, P::Number) where {uType, tType}
    throw(ErrorException("Lindemann rate constants not implemented yet."))
end

function has_conditions(::PrecalculatedLindemannCalculator, symbols::Vector{Symbol})
    return all([sym in [:T, :P] for sym in symbols])
end

allows_continuous(::PrecalculatedLindemannCalculator) = true