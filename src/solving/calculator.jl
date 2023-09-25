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
struct DummyKineticCalculator{kmType, uType, tType} <: AbstractKineticCalculator
    rates::Vector{uType}
    k_max::kmType
    t_unit::String
    t_mult::tType
end

"""
    calculator = DummyKineticCalculator(rd, rates[, k_max, t_unit])

Outer constructor method for placeholder kinetic calculator.
"""
function DummyKineticCalculator(rd::RxData, rates::Vector{uType}; 
        k_max::Union{Nothing, uType}=nothing, t_unit::String="s") where {uType <: AbstractFloat}
    if length(rates) != rd.nr
        throw(ArgumentError("Number of rates ($(length(rates))) does not match number of reactions in `RxData` ($(rd.nr))"))
    end
    t_mult = tconvert(t_unit, "s")
    return DummyKineticCalculator(rates, k_max, t_unit, t_mult)
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
struct PrecalculatedArrheniusCalculator{kmType, uType, tType} <: AbstractKineticCalculator
    Ea::Vector{uType}
    A::Vector{uType}
    k_max::kmType
    t_unit::String
    t_mult::tType
end

"""
    calculator = PrecalculatedArrheniusCalculator(rd, rates[, k_max, t_unit])

Outer constructor method for Arrhenius theory kinetic calculator.
"""
function PrecalculatedArrheniusCalculator(rd::RxData, Ea::Vector{uType}, A::Vector{uType}; 
        k_max::Union{Nothing, uType}=nothing, t_unit::String="s") where {uType <: AbstractFloat}
    if length(Ea) != rd.nr || length(A) != rd.nr
        throw(ArgumentError("Number of parameters (Ea: $(length(Ea)), A: $(length(A))) does not match number of reactions in `RxData` ($(rd.nr))"))
    end
    t_mult = tconvert(t_unit, "s")
    return PrecalculatedArrheniusCalculator(Ea, A, k_max, t_unit, t_mult)
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
    R = 8.314462618 # Gas constant (J/K/mol)
    N_A = 6.02214076e23 # Avogadro constant (/mol)

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
struct PrecalculatedLindemannCalculator{kmType, uType, tType} <: AbstractKineticCalculator
    Ea::Vector{uType}
    A_0::Vector{uType}
    A_inf::Vector{uType}
    k_max::kmType
    t_unit::String
    t_mult::tType
end

"""
    calculator = PrecalculatedLindemannCalculator(rd, rates[, k_max, t_unit])

Outer constructor method for pressure-dependent Arrhenius theory kinetic calculator.
"""
function PrecalculatedLindemannCalculator(rd::RxData, Ea::Vector{uType}, A::Vector{uType}; 
        k_max::Union{Nothing, uType}=nothing, t_unit::String="s") where {uType <: AbstractFloat}
    if length(rates) != rd.nr
        throw(ArgumentError("Number of rates ($(length(rates))) does not match number of reactions in `RxData` ($(rd.nr))"))
    end
    t_mult = tconvert(t_unit, "s")
    return PrecalculatedLindemannCalculator(Ea, A, k_max, t_unit, t_mult)
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