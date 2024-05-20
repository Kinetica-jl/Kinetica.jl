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
    solve_variable_condition!(profile<:AbstractGradientProfile, pars::ODESimulationParams[, sym=nothing, reset=false, solver, solve_kwargs])

Generates a solution for the specified gradient-variable condition profile.

For gradient-based profiles, this requires constructing an
ODEProblem around their MTK-derived symbolic gradient expressions
and solving over the timespan in `pars`.

The ODE solver and the arguments passed to the `solve()` call
can be controlled with the `solver` and `solve_kwargs` arguments
respectively. If `sym` is passed a `Symbol`, this will bind the
solution result to that symbol in the underlying `ODESolution`.
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
        @named profile_sys = ODESystem([D(X_sym) ~ profile.grad(t, profile)], t)

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


# -----------------------------------------------
# NullGradientProfile definition
# -----------------------------------------------
mutable struct NullGradientProfile{uType, tType} <: AbstractGradientProfile
    grad::Function
    X_start::uType
    t_end::tType
    tstops::Vector{tType}
    sol
end

"""
    NullGradientProfile(; X_start, t_end)

Container for null gradient profile data and gradient function.

This condition profile should only be used for debugging,
as it has a condition gradient function which always
returns zero (i.e. there is no condition change from
`X_start`). If only this constant condition is required, 
`StaticODESolve` should always be used with a
`StaticConditionProfile` instead of a `VariableODESolve`
with this condition profile.

Contains fields for:
* Condition gradient function (`grad`)
* Initial value of condition (`X_start`)
* Time to stop calculation (`t_end`)
* Times for the ODE solver to ensure calculation at (`tstops`)
* Profile solution, constructed by call to `solve_variable_condition!` (`sol`)
"""
function NullGradientProfile(;
    X_start::uType,
    t_end::tType,
) where {uType <: AbstractFloat, tType <: AbstractFloat}

    tstops = [t_end]
    return NullGradientProfile(_grad_NullGradientProfile, X_start, t_end, tstops, nothing)
end

function _grad_NullGradientProfile(t, profile::NullGradientProfile)
    return typeof(profile.X_start)(0.0)
end

function create_discrete_tstops!(profile::NullGradientProfile, ts_update::AbstractFloat)
    if ts_update > profile.t_end throw(ArgumentError("Error defining tstops, `ts_update` is too large.")) end
    profile.tstops = collect(0.0:ts_update:profile.t_end)
end


# -----------------------------------------------
# LinearGradientProfile definition
# -----------------------------------------------
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
    LinearGradientProfile(; rate, X_start, X_end)

Container for linear condition ramp profile data and condition gradient function.

This condition profile represents a linear condition
increase/decrease from `X_start` to `X_end`. Determines 
the simulation end time from the provided conditions
and gradient, then constructs the condition gradient 
function (which returns `rate` for every timestep).

Contains fields for:
* Condition gradient function (`grad`)
* Rate of change of condition (`rate`)
* Initial value of condition (`X_start`)
* Final value of condition (`X_end`)
* Time to stop calculation (`t_end`)
* Times for the ODE solver to ensure calculation at (`tstops`)
* Profile solution, constructed by call to `solve_variable_condition!` (`sol`)
"""
function LinearGradientProfile(;
    rate::uType,
    X_start::uType,
    X_end::uType
) where {uType <: AbstractFloat}

    if (X_end < X_start && rate > 0) || (X_end > X_start && rate < 0)
        error("Impossible condition ramp defined. Check heating rates have the correct signs.")
    end
    t_end = (X_end - X_start)/rate
    tstops = [t_end]

    return LinearGradientProfile(_grad_LinearGradientProfile, rate, X_start, X_end, t_end, tstops, nothing)
end

function _grad_LinearGradientProfile(t, profile::LinearGradientProfile)
    return typeof(profile.X_start)(
        ((t <= profile.t_end) * profile.rate) +
        ((t > profile.t_end) * 0.0)
    )
end

function create_discrete_tstops!(profile::LinearGradientProfile, ts_update::AbstractFloat)
    if ts_update > profile.t_end throw(ArgumentError("Error defining tstops, `ts_update` is too large.")) end
    profile.tstops = create_savepoints(0.0, profile.t_end, ts_update)
end


# -----------------------------------------------
# DoubleRampGradientProfile definition
# -----------------------------------------------
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
    t_startr1::tType
    t_endr1::tType
    t_startr2::tType
    t_endr2::tType
    t_blend::tType
    t_end::tType
    tstops::Vector{tType}
    sol
end

"""
    DoubleRampGradientProfile(; X_start, t_start_plateau, rate1, X_mid, t_mid_plateau, rate2, X_end, t_end_plateau[, t_blend])

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
        error("Impossible condition ramp defined. Check heating rates have the correct signs.")
    end

    t_startr1 = t_start_plateau
    t_endr1 = tType(t_startr1+ ((X_mid - X_start)/rate1))
    t_startr2 = t_endr1 + t_mid_plateau
    t_endr2 = tType(t_startr2 + ((X_end - X_mid)/rate2))
    t_end = t_endr2 + t_end_plateau

    if isnothing(t_blend)
        tstops = [t_startr1, t_endr1, t_startr2, t_endr2, t_end]
        return DoubleRampGradientProfile(
            _grad_DoubleRampGradientProfile,
            rate1, rate2, X_start, X_mid, X_end, 
            t_start_plateau, t_mid_plateau, t_end_plateau, 
            t_startr1, t_endr1, t_startr2, t_endr2, tType(0.0), 
            t_end, tstops, nothing)
    else
        tstops = [
            t_startr1-t_blend, t_startr1+t_blend,
            t_endr1-t_blend, t_endr1+t_blend,
            t_startr2-t_blend, t_startr2+t_blend,
            t_endr2-t_blend, t_endr2+t_blend,
            t_end
        ]
        return DoubleRampGradientProfile(
            _grad_DoubleRampGradientProfile_blended,
            rate1, rate2, X_start, X_mid, X_end,
            t_start_plateau, t_mid_plateau, t_end_plateau,
            t_startr1, t_endr1, t_startr2, t_endr2, t_blend,
            t_end, tstops, nothing)
    end
end

function _grad_DoubleRampGradientProfile(t, p::DoubleRampGradientProfile)
    return typeof(p.X_start)(
        ((t < p.t_startr1) * 0.0) +
        ((t >= p.t_startr1 && t < p.t_endr1) * p.rate1) +
        ((t >= p.t_endr1 && t < p.t_startr2) * 0.0) +
        ((t >= p.t_startr2 && t < p.t_endr2) * p.rate2) +
        ((t >= p.t_endr2) * 0.0) 
    )
end

function _grad_DoubleRampGradientProfile_blended(t, p::DoubleRampGradientProfile)
    return typeof(p.X_start)(
        ((t < p.t_startr1-p.t_blend) * 0.0) +
        ((t >= p.t_startr1-p.t_blend && t < p.t_startr1+p.t_blend) * (p.rate1*(t-p.t_startr1-p.t_blend)/(2*p.t_blend) + p.rate1)) +
        ((t >= p.t_startr1+p.t_blend && t < p.t_endr1-p.t_blend) * p.rate1) +
        ((t >= p.t_endr1-p.t_blend && t < p.t_endr1+p.t_blend) * (-p.rate1*(t-p.t_endr1-p.t_blend)/(2*p.t_blend))) +
        ((t >= p.t_endr1+p.t_blend && t < p.t_startr2-p.t_blend) * 0.0) +
        ((t >= p.t_startr2-p.t_blend && t < p.t_startr2+p.t_blend) * (p.rate2*(t-p.t_startr2-p.t_blend)/(2*p.t_blend) + p.rate2)) +
        ((t >= p.t_startr2+p.t_blend && t < p.t_endr2-p.t_blend) * p.rate2) +
        ((t >= p.t_endr2-p.t_blend && t < p.t_endr2+p.t_blend) * (-p.rate2*(t-p.t_endr2-p.t_blend)/(2*p.t_blend))) +
        ((t >= p.t_endr2+p.t_blend) * 0.0)
    )
end

function create_discrete_tstops!(profile::DoubleRampGradientProfile, ts_update::AbstractFloat)
    if ts_update > profile.t_end throw(ArgumentError("Error defining tstops, `ts_update` is too large.")) end

    profile.tstops = reduce(vcat, [
        [0.0],
        create_savepoints(profile.t_startr1-profile.t_blend, profile.t_endr1+profile.t_blend, ts_update),
        create_savepoints(profile.t_startr2-profile.t_blend, profile.t_endr2+profile.t_blend, ts_update),
        [profile.t_end]
    ])
end