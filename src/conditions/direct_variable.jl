"""
Definitions for all directly-variable condition profile structs and their outer constructor methods.

Directly-variable conditions are those that are evolved
in time through a direct function, `f(t)`. These profiles
do not need to be solved by an ODE solver and have an
exact functional form for their condition at any given time.

Usable with any condition, which are referred to by the
name `X`.

For compatibility, all direct profile structs must implement the following fields:

* `f<:Function`
* `X_start<:AbstractFloat`
* `t_end<:AbstractFloat`
* `tstops::Vector{<:AbstractFloat}`
* `sol`
"""


"""
    solve_variable_condition!(profile<:AbstractDirectProfile, pars[, reset, solve_kwargs...])

Generates a solution for the specified directly-variable condition profile.

For profiles with direct functions, this requires calculating values
for the specified `pars.tspan` and wrapping them within an `ODESolution`
for compatibility with other interfaces (plotting, interpolation, etc.).
"""
function solve_variable_condition!(profile::pType, pars::ODESimulationParams;
    reset=false, solve_kwargs...) where {pType <: AbstractDirectProfile}
    if isnothing(profile.sol) || reset            
        save_interval = isnothing(pars.save_interval) ? pars.tspan[2]/1000 : pars.save_interval
        t = create_savepoints(pars.tspan[1], pars.tspan[2], save_interval)
        u = [[profile.f(tp)] for tp in t]
        profile.sol = ODESolution{typeof(profile.X_start), 1}(
            u,
            nothing,
            nothing,
            t,
            nothing,
            DummyODEProblem(),
            nothing,
            SciMLBase.LinearInterpolation(t, u),
            false,
            0,
            nothing,
            nothing,
            ReturnCode.Default
        )
    end
    return
end


"""
Container for null direct profile data and condition function.

This condition profile should only be used for debugging,
as it has a condition function which always returns the
initial condition. If only this constant condition is 
required, the regular `ODESolve` should always be used 
instead of an `ODEConditionSolve` with this condition 
profile.

Contains fields for:
* Condition function (`f`)
* Initial value of condition (`X_start`)
* Time to stop calculation (`t_end`)
* Times for the ODE solver to ensure calculation at (`tstops`)
* Profile solution, constructed by call to `solve_variable_condition!` (`sol`)
"""
mutable struct NullDirectProfile{uType, tType} <: AbstractDirectProfile
    f::Function
    X_start::uType
    t_end::tType
    tstops::Vector{tType}
    sol
end

"""
    condition_profile = NullDirectProfile(; X, t_end)

Outer constructor for null condition direct profile.

Should only be used for testing purposes (see struct
documentation).
"""
function NullDirectProfile(;
    X::uType,
    t_end::tType,
) where {uType <: AbstractFloat, tType <: AbstractFloat}
    function f(t)
        return X
    end

    tstops = [t_end]

    return NullDirectProfile(f, X, t_end, tstops, nothing)
end

function create_discrete_tstops(profile::NullDirectProfile, ts_update::AbstractFloat)
    if ts_update > profile.t_end throw(ArgumentError("Error defining tstops, `ts_update` is too large.")) end
    profile.tstops = collect(0.0:ts_update:profile.t_end)
end


"""
Container for linear condition ramp profile data and condition function.

This condition profile represents a linear condition
increase/decrease from `X_start` to `X_end`.

Contains fields for:
* Condition function (`f`)
* Rate of change of condition (`rate`)
* Initial value of condition (`X_start`)
* Final value of condition (`X_end`)
* Time to stop calculation (`t_end`)
* Times for the ODE solver to ensure calculation at (`tstops`)
* Profile solution, constructed by call to `solve_variable_condition!` (`sol`)
"""
mutable struct LinearDirectProfile{uType, tType} <: AbstractDirectProfile
    f::Function
    rate::uType
    X_start::uType
    X_end::uType
    t_end::tType
    tstops::Vector{tType}
    sol
end

"""
    condition_profile = LinearDirectProfile(; rate, X_start, X_end)

Outer constructor for linear condition ramp direct profile.

Determines the simulation end time from the provided conditions
and rate, then constructs the condition function (which is a linear
y = mx + c function).
"""
function LinearDirectProfile(;
    rate::uType,
    X_start::uType,
    X_end::uType
) where {uType <: AbstractFloat}
    t_end = (X_end - X_start)/rate

    function f(t)
        return uType(
            ((t <= 0.0) * X_start) +
            ((t > 0.0 && t <= t_end) * (X_start + (rate * t))) + 
            ((t > t_end) * X_end)
        )
    end

    tstops = [t_end]

    return LinearDirectProfile(f, rate, X_start, X_end, t_end, tstops, nothing)
end

function create_discrete_tstops(profile::LinearDirectProfile, ts_update::AbstractFloat)
    if ts_update > profile.t_end throw(ArgumentError("Error defining tstops, `ts_update` is too large.")) end
    profile.tstops = create_savepoints(0.0, profile.t_end, ts_update)
end