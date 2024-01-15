"""
Definitions for all gradient-variable condition profile structs and their outer constructor methods.

Gradient-variable conditions are those that are evolved
in time through a gradient function. These must be solved
e.g. with an ODE solver to elucidate the whole profile.

Usable with any condition, which are referred to by the
name `X`.

For compatibility, all gradient profile structs must implement the following fields:

* `grad<:Function`
* `X_start<:AbstractFloat`
* `t_end<:AbstractFloat`
* `tstops::Vector{<:AbstractFloat}`
* `sol`
"""


"""
    solve_variable_condition!(profile<:AbstractGradientProfile, pars[, reset, solver, solve_kwargs])

Generates a solution for the specified gradient-variable condition profile.

For gradient-based profiles, this requires constructing an
ODEProblem around their MTK-derived symbolic gradient expressions
and solving over the timespan in `pars`.

The ODE solver and the arguments passed to the `solve()` call
can be controlled with the `solver` and `solve_kwargs` arguments
respectively.
"""
function solve_variable_condition!(profile::pType, pars::ODESimulationParams;
    sym=nothing, reset=false, solver, solve_kwargs) where {pType <: AbstractGradientProfile}
    if isnothing(profile.sol) || reset
        @variables t 
        if isnothing(sym) 
            X_sym = only(@variables(X(t)))
        else
            X_sym = only(@variables($sym(t)))
        end
        D = Differential(t)
        @named profile_sys = ODESystem([D(X_sym) ~ profile.grad(t)], t)

        u0map = [Pair(X_sym, profile.X_start)]
        profile_solve_kwargs = deepcopy(solve_kwargs)
        if :tstops in keys(profile_solve_kwargs)
            profile_solve_kwargs[:tstops] = sort(vcat(profile_solve_kwargs[:tstops], profile.tstops))
        else
            profile_solve_kwargs[:tstops] = profile.tstops
        end
        if !(:saveat in keys(profile_solve_kwargs))
            save_interval = isnothing(pars.save_interval) ? pars.tspan[2]/1000 : pars.save_interval
            profile_solve_kwargs[:saveat] = sort(vcat(create_savepoints(pars.tspan[1], pars.tspan[2], save_interval), profile.tstops))
        end

        prob = ODEProblem(profile_sys, u0map, pars.tspan)
        profile.sol = solve(prob, solver; profile_solve_kwargs...)
    end
    return
end


"""
Container for null gradient profile data and gradient function.

This condition profile should only be used for debugging,
as it has a condition gradient function which always
returns zero (i.e. there is no condition change from
`X_start`). If only this constant condition is required, the
regular `ODESolve` should always be used instead of an
`ODEConditionSolve` with this condition profile.

Contains fields for:
* Condition gradient function (`grad`)
* Initial value of condition (`X_start`)
* Time to stop calculation (`t_end`)
* Times for the ODE solver to ensure calculation at (`tstops`)
* Profile solution, constructed by call to `solve_variable_condition!` (`sol`)
"""
mutable struct NullGradientProfile{uType, tType} <: AbstractGradientProfile
    grad::Function
    X_start::uType
    t_end::tType
    tstops::Vector{tType}
    sol
end

"""
    condition_profile = NullGradientProfile(; X, t_end)

Outer constructor for null condition gradient profile.

Should only be used for testing purposes (see struct
documentation).
"""
function NullGradientProfile(;
    X::uType,
    t_end::tType,
) where {uType <: AbstractFloat, tType <: AbstractFloat}
    function grad(t)
        return 0.0
    end

    tstops = [t_end]

    return NullGradientProfile(grad, X, t_end, tstops, nothing)
end

function create_discrete_tstops(profile::NullGradientProfile, ts_update::AbstractFloat)
    if ts_update > profile.t_end throw(ArgumentError("Error defining tstops, `ts_update` is too large.")) end
    profile.tstops = collect(0.0:ts_update:profile.t_end)
end



"""
Container for linear condition ramp profile data and condition gradient function.

This condition profile represents a linear condition
increase/decrease from `X_start` to `X_end`.

Contains fields for:
* Condition gradient function (`grad`)
* Rate of change of condition (`rate`)
* Initial value of condition (`X_start`)
* Final value of condition (`X_end`)
* Time to stop calculation (`t_end`)
* Times for the ODE solver to ensure calculation at (`tstops`)
* Profile solution, constructed by call to `solve_variable_condition!` (`sol`)
"""
mutable struct LinearGradientProfile{uType, tType} <: AbstractGradientProfile
    grad::Function
    rate::uType
    X_start::uType
    X_end::uType
    t_end::tType
    tstops::Vector{tType}
    sol
end

"""
    condition_profile = LinearGradientProfile(; rate, X_start, X_end)

Outer constructor for linear condition ramp gradient profile.

Determines the simulation end time from the provided conditions
and gradient, then constructs the condition gradient function 
(which returns `rate` for every timestep).
"""
function LinearGradientProfile(;
    rate::uType,
    X_start::uType,
    X_end::uType
) where {uType <: AbstractFloat}

    if (X_end < X_start && rate > 0) || (X_end > X_start && rate < 0)
        error("Impossible temperature ramp defined. Check heating rates have the correct signs.")
    end

    t_end = (X_end - X_start)/rate

    function grad(t)
        return ((t <= t_end) * rate) +
               ((t > t_end) * 0.0)
    end

    tstops = [t_end]

    return LinearGradientProfile(grad, rate, X_start, X_end, t_end, tstops, nothing)
end

function create_discrete_tstops(profile::LinearGradientProfile, ts_update::AbstractFloat)
    if ts_update > profile.t_end throw(ArgumentError("Error defining tstops, `ts_update` is too large.")) end
    profile.tstops = create_savepoints(0.0, profile.t_end, ts_update)
end


"""
Container for double condition ramp profile data and condition gradient function.

This condition profile represents two condition ramps with
adjustable condition plateaus before, after and in between
the ramps, i.e.

                  ------   X_mid
          rate1  /      \\
                /        \\  rate2
    X_start ----          \\
                           ----- X_end

The profile starts at `X_start` and maintains that value for
`t_start_plateau`. The condition then ramps with gradient
`rate1` to condition value `X_mid`. This value is maintained
for `t_mid_plateau`. The condition then ramps with gradient
`rate2` to condition value `X_end`. This value is maintained
for `t_end_plateau` until the calculated time `t_end`.

To smooth out gradient discontinuities, a blending time `t_blend`
can be passed to linearly interploate between plateaus and
ramps, forming a smooth function of time. Larger values of
`t_blend` yield smoother functions but decrease accuracy of
the ramps, so should be used carefully.
"""
mutable struct DoubleRampGradientProfile{uType, tType} <: AbstractGradientProfile
    grad::Function
    rate1::uType
    rate2::uType
    X_start::uType
    X_mid::uType
    X_end::uType
    t_start_plateau::tType
    t_mid_plateau::tType
    t_end_plateau::tType
    t_blend::tType
    t_end::tType
    tstops::Vector{tType}
    sol
end

"""
    condition_profile = DoubleRampGradientProfile(; X_start, t_start_plateau, rate1, X_mid, t_mid_plateau, rate2, X_end, t_end_plateau[, t_blend])

Outer constructor for double condition ramp gradient profile.

Determines the simulation end time from the provided conditions
and gradients, then constructs the condition gradient function 
(which returns `rate` for every timestep).

If `t_blend` is passed, constructs a smooth approximation of
the otherwise discontinuous gradient function.
"""
function DoubleRampGradientProfile(;
    X_start::uType,
    t_start_plateau::tType,
    rate1::uType,
    X_mid::uType,
    t_mid_plateau::tType,
    rate2::uType,
    X_end::uType,
    t_end_plateau::tType,
    t_blend::Union{Nothing, tType} = nothing
) where {uType <: AbstractFloat, tType <: AbstractFloat}
    
    if (X_mid > X_start && rate1 < 0) || (X_mid < X_start && rate1 > 0) ||
        (X_end > X_mid && rate2 < 0) || (X_end < X_mid && rate2 > 0)
        error("Impossible temperature ramp defined. Check heating rates have the correct signs.")
    end

    t_startr1 = t_start_plateau
    t_endr1 = tType(t_startr1+ ((X_mid - X_start)/rate1))
    t_startr2 = t_endr1 + t_mid_plateau
    t_endr2 = tType(t_startr2 + ((X_end - X_mid)/rate2))
    t_end = t_endr2 + t_end_plateau

    function grad(t)
        return uType(
            ((t < t_startr1) * 0.0) +
            ((t >= t_startr1 && t < t_endr1) * rate1) +
            ((t >= t_endr1 && t < t_startr2) * 0.0) +
            ((t >= t_startr2 && t < t_endr2) * rate2) +
            ((t >= t_endr2) * 0.0) 
        )
    end

    function grad_blend(t)
        return uType(
            ((t < t_startr1-t_blend) * 0.0) +
            ((t >= t_startr1-t_blend && t < t_startr1+t_blend) * (rate1*(t-t_startr1-t_blend)/(2*t_blend) + rate1)) +
            ((t >= t_startr1+t_blend && t < t_endr1-t_blend) * rate1) +
            ((t >= t_endr1-t_blend && t < t_endr1+t_blend) * (-rate1*(t-t_endr1-t_blend)/(2*t_blend))) +
            ((t >= t_endr1+t_blend && t < t_startr2-t_blend) * 0.0) +
            ((t >= t_startr2-t_blend && t < t_startr2+t_blend) * (rate2*(t-t_startr2-t_blend)/(2*t_blend) + rate2)) +
            ((t >= t_startr2+t_blend && t < t_endr2-t_blend) * rate2) +
            ((t >= t_endr2-t_blend && t < t_endr2+t_blend) * (-rate2*(t-t_endr2-t_blend)/(2*t_blend))) +
            ((t >= t_endr2+t_blend) * 0.0)
        )
    end

    if isnothing(t_blend)
        tstops = [t_startr1, t_endr1, t_startr2, t_endr2, t_end]
        return DoubleRampGradientProfile(
            grad, rate1, rate2, X_start, X_mid, X_end, 
            t_start_plateau, t_mid_plateau, t_end_plateau, 
            tType(0.0), t_end, tstops, nothing)
    else
        tstops = [
            t_startr1-t_blend, t_startr1+t_blend,
            t_endr1-t_blend, t_endr1+t_blend,
            t_startr2-t_blend, t_startr2+t_blend,
            t_endr2-t_blend, t_endr2+t_blend,
            t_end
        ]
        return DoubleRampGradientProfile(
            grad_blend, rate1, rate2, X_start, X_mid, X_end,
            t_start_plateau, t_mid_plateau, t_end_plateau,
            t_blend, t_end, tstops, nothing)
    end
end

function create_discrete_tstops(profile::DoubleRampGradientProfile, ts_update::AbstractFloat)
    if ts_update > profile.t_end throw(ArgumentError("Error defining tstops, `ts_update` is too large.")) end

    tType = eltype(profile.tstops)
    t_startr1 = profile.t_start_plateau
    t_endr1 = tType(t_startr1 + ((profile.X_mid - profile.X_start)/profile.rate1))
    t_startr2 = t_endr1 + profile.t_mid_plateau
    t_endr2 = tType(t_startr2 + ((profile.X_end - profile.X_mid)/profile.rate2))
    profile.tstops = reduce(vcat, [
        [0.0],
        create_savepoints(t_startr1-profile.t_blend, t_endr1+profile.t_blend, ts_update),
        create_savepoints(t_startr2-profile.t_blend, t_endr2+profile.t_blend, ts_update),
        [profile.t_end]
    ])
end