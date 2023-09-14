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
    solve_variable_condition!(profile<:AbstractGradientProfile, pars[, reset, solve_kwargs...])

Generates a solution for the specified gradient-variable condition profile.

For gradient-based profiles, this requires constructing an
ODEProblem around their MTK-derived symbolic gradient expressions
and solving over the timespan in `pars`.
"""
function solve_variable_condition!(profile::pType, pars::ODESimulationParams;
    reset=false, solve_kwargs...) where {pType <: AbstractGradientProfile}
    if isnothing(profile.sol) || reset
        @variables t X(t)
        D = Differential(t)
        @named profile_sys = ODESystem([D(X) ~ profile.grad(t)], t)
        u0map = [Pair(X, profile.X_start)]
        prob = ODEProblem(profile_sys, u0map, pars.tspan)
        profile.sol = solve(prob; solve_kwargs...)
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


# Old temperature ramp code, needs to be converted for general conditions

# """
# Container for simple double temperature ramp profile data and temperature gradient function.

# This temperature profile steps through an initial plateau, 
# followed by a discontinuous ramp to a middle temperature,
# followed by another plateau, followed by another discontinuous
# ramp to a final temperature, ending after a final plateau.

# Note that while both `rate1` and `rate2` are defined as heating
# rates, a negative value for either/both indicates a downward
# temperature ramp, i.e. a cooling.

# Contains fields for:
# * Temperature gradient function (`Tgrad`)
# * 1st heating rate (`rate1`)
# * 2nd heating rate (`rate2`)
# * Initial temperature (`T_start`)
# * Middle plateau temperature (`T_mid`)
# * Final temperature (`T_end`)
# * Maximum temperature (`T_max`)
# * Length of time to spend at initial temperature plateau (`t_startplat`)
# * Length of time to spend at maximum temperature plateau (`t_midplat`)
# * Length of time to spend at final temperature plateau (`t_endplat`)
# * Time to start 1st heating ramp (`t_startr1`)
# * Time to stop 1st heating ramp (`t_endr1`)
# * Time to start 2nd heating ramp (`t_startr2`)
# * Time to stop 2nd heating ramp (`t_endr2`)
# * Time to stop calculation (`t_end`)
# * Whether discrete rate updates should be enabled (`discrete_rates`)
# * Timestep for discrete rate updates (`discrete_rate_timestep`)
# * Times for the ODE solver to ensure calculation at (`tstops`)
# """
# struct SimpleDoubleRampTprofile{uType, tType} <: AbstractTprofile
#     Tgrad::Function
#     rate1::uType
#     rate2::uType
#     T_start::uType
#     T_mid::uType
#     T_end::uType
#     T_max::uType
#     t_startplat::tType
#     t_midplat::tType
#     t_endplat::tType
#     t_startr1::tType
#     t_endr1::tType
#     t_startr2::tType
#     t_endr2::tType
#     t_end::tType
#     discrete_rates::Bool
#     discrete_rate_timestep::Union{tType, Nothing}
#     tstops::Vector{tType}
# end

# """
#     Tprofile = SimpleDoubleRampTprofile(; rate1=rate1, rate2=rate2, 
#         T_start=T_start, T_mid=T_mid, T_end=T_end, t_startplat=t_startplat,
#         t_midplat=t_midplat, t_endplat=t_endplat[, discrete_rate_timestep=discrete_rate_timestep])

# Outer constructor for simple double temperature ramp temperature profile.

# Creates all of the required absolute times from the provided 
# relative times, then constructs the temperature gradient function
# for the specified simple temperature ramp profile.

# Also creates the `tstops` array for passing to adaptive timestepping
# ODE solvers.
# """
# function SimpleDoubleRampTprofile(;
#     rate1::uType,
#     rate2::uType,
#     T_start::uType,
#     T_mid::uType,
#     T_end::uType,
#     t_startplat::tType,
#     t_midplat::tType,
#     t_endplat::tType,
#     discrete_rate_timestep::Union{tType, Nothing} = nothing
# ) where {uType <: AbstractFloat, tType <: AbstractFloat}
#     if (T_mid < T_start && rate1 > 0) || (T_mid > T_end && rate2 < 0)
#         error("Impossible temperature ramps defined. Check heating rates have the correct signs.")
#     end

#     t_startr1 = t_startplat
#     t_endr1 = tType(t_startr1+ ((T_mid - T_start)/rate1))
#     t_startr2 = t_endr1 + t_midplat
#     t_endr2 = tType(t_startr2 + ((T_end - T_mid)/rate2))
#     t_end = t_endr2 + t_endplat
#     T_max = maximum([T_start, T_mid, T_end])

#     function Tgrad(t)
#         return uType(
#             ((t < t_startr1) * 0.0) +
#             ((t >= t_startr1 && t < t_endr1) * rate1) +
#             ((t >= t_endr1 && t < t_startr2) * 0.0) +
#             ((t >= t_startr2 && t < t_endr2) * rate2) +
#             ((t >= t_endr2) * 0.0) 
#         )
#     end

#     discrete = false
#     if !isnothing(discrete_rate_timestep)
#         if discrete_rate_timestep > t_end throw(ArgumentError("Error defining tstops, discrete_rate_timestep is too large.")) end
#         tstops = reduce(vcat,
#             [
#                 [tType(0.0)],
#                 create_savepoints(t_startr1, t_endr1, discrete_rate_timestep),
#                 create_savepoints(t_startr2, t_endr2, discrete_rate_timestep),
#                 [t_end]
#             ])
#         discrete = true
#     else
#         tstops = [
#             t_startr1, t_endr1,
#             t_startr2, t_endr2,
#             t_end
#         ]
#     end

#     return SimpleDoubleRampTprofile(
#         Tgrad, rate1, rate2, T_start, T_mid, T_end, T_max, t_startplat, t_midplat,
#         t_endplat, t_startr1, t_endr1, t_startr2, t_endr2, t_end, discrete,
#         discrete_rate_timestep, tstops
#     )
# end


# """
# Container for smooth double temperature ramp profile data and temperature gradient function.

# This temperature profile steps through an initial plateau, 
# followed by a smooth ramp up to a maximum temperature,
# followed by another plateau, followed by a smooth ramp down
# to a final temperature, ending after a final plateau.

# Note that while both `rate1` and `rate2` are defined as heating
# rates, a negative value for either/both indicates a downward
# temperature ramp, i.e. a cooling.

# Contains fields for:
# * Temperature gradient function (`Tgrad`)
# * 1st heating rate (`rate1`)
# * 2nd heating rate (`rate2`)
# * Initial temperature (`T_start`)
# * Middle plateau temperature (`T_mid`)
# * Final temperature (`T_end`)
# * Maximum temperature (`T_max`)
# * Length of time to spend at initial temperature plateau (`t_startplat`)
# * Length of time to spend at maximum temperature plateau (`t_midplat`)
# * Length of time to spend at final temperature plateau (`t_endplat`)
# * Blending time to be taken out of each plateau (`t_blend`)
# * Time to start 1st heating ramp (`t_startr1`)
# * Time to stop 1st heating ramp (`t_endr1`)
# * Time to start 2nd heating ramp (`t_startr2`)
# * Time to stop 2nd heating ramp (`t_endr2`)
# * Time to stop calculation (`t_end`)
# * Whether discrete rate updates should be enabled (`discrete_rates`)
# * Timestep for discrete rate updates (`discrete_rate_timestep`)
# * Times for the ODE solver to ensure calculation at (`tstops`)
# """
# struct SmoothDoubleRampTprofile{uType,tType} <: AbstractTprofile
#     Tgrad::Function
#     rate1::uType
#     rate2::uType
#     T_start::uType
#     T_mid::uType
#     T_end::uType
#     T_max::uType
#     t_startplat::tType
#     t_midplat::tType
#     t_endplat::tType
#     t_blend::tType
#     t_startr1::tType
#     t_endr1::tType
#     t_startr2::tType
#     t_endr2::tType
#     t_end::tType
#     discrete_rates::Bool
#     discrete_rate_timestep::Union{tType, Nothing}
#     tstops::Vector{tType}
# end

# """
#     Tprofile = SmoothDoubleRampTprofile(; heatrate=heatrate, coolrate=coolrate, 
#         T_start=T_start, T_max=T_max, T_end=T_end, t_startplateau=t_startplateau,
#         t_midplateau=t_midplateau, t_endplateau=t_endplateau, t_blend=t_blend
#         [, discrete_rate_timestep=discrete_rate_timestep])

# Outer constructor for smooth double temperature ramp temperature profile.

# Creates all of the required absolute times from the provided 
# relative times, then constructs the temperature gradient function
# for the specified smooth temperature ramp profile.

# Also creates the `tstops` array for passing to adaptive timestepping
# ODE solvers.
# """
# function SmoothDoubleRampTprofile(;
#     rate1::uType,
#     rate2::uType,
#     T_start::uType,
#     T_mid::uType,
#     T_end::uType,
#     t_startplat::tType,
#     t_midplat::tType,
#     t_endplat::tType,
#     t_blend::tType,
#     discrete_rate_timestep::Union{tType, Nothing} = nothing
# ) where {uType <: AbstractFloat, tType <: AbstractFloat}
#     if (T_mid < T_start && rate1 > 0) || (T_mid > T_end && rate2 < 0)
#         error("Impossible temperature ramps defined. Check heating rates have the correct signs.")
#     end
#     if t_blend < 0.0
#         error("t_blend must be a positive value.")
#     end

#     t_startr1 = t_startplat
#     t_endr1 = tType(t_startr1+ ((T_mid - T_start)/rate1))
#     t_startr2 = t_endr1 + t_midplat
#     t_endr2 = tType(t_startr2 + ((T_end - T_mid)/rate2))
#     t_end = t_endr2 + t_endplat
#     T_max = maximum([T_start, T_mid, T_end])

#     function Tgrad(t)
#         return uType(
#             ((t < t_startr1-t_blend) * 0.0) +
#             ((t >= t_startr1-t_blend && t < t_startr1+t_blend) * (rate1*(t-t_startr1-t_blend)/(2*t_blend) + rate1)) +
#             ((t >= t_startr1+t_blend && t < t_endr1-t_blend) * rate1) +
#             ((t >= t_endr1-t_blend && t < t_endr1+t_blend) * (-rate1*(t-t_endr1-t_blend)/(2*t_blend))) +
#             ((t >= t_endr1+t_blend && t < t_startr2-t_blend) * 0.0) +
#             ((t >= t_startr2-t_blend && t < t_startr2+t_blend) * (rate2*(t-t_startr2-t_blend)/(2*t_blend) + rate2)) +
#             ((t >= t_startr2+t_blend && t < t_endr2-t_blend) * rate2) +
#             ((t >= t_endr2-t_blend && t < t_endr2+t_blend) * (-rate2*(t-t_endr2-t_blend)/(2*t_blend))) +
#             ((t >= t_endr2+t_blend) * 0.0)
#         )
#     end

#     discrete = false
#     if !isnothing(discrete_rate_timestep)
#         if discrete_rate_timestep > t_end throw(ArgumentError("Error defining tstops, discrete_rate_timestep is too large.")) end
#         tstops = reduce(vcat,
#             [
#                 [tType(0.0)],
#                 create_savepoints(t_startr1-t_blend, t_endr1+t_blend, discrete_rate_timestep),
#                 create_savepoints(t_startr2-t_blend, t_endr2+t_blend, discrete_rate_timestep),
#                 [t_end]
#             ])
#         discrete = true
#     else
#         tstops = [
#             t_startr1-t_blend, t_startr1, t_startr1+t_blend, 
#             t_endr1-t_blend, t_endr1, t_endr1+t_blend,
#             t_startr2-t_blend, t_startr2, t_startr2+t_blend,
#             t_endr2-t_blend, t_endr2, t_endr2+t_blend
#         ]
#     end

#     return SmoothDoubleRampTprofile(
#         Tgrad, rate1, rate2, T_start, T_mid, T_end, T_max, t_startplat, t_midplat,
#         t_endplat, t_blend, t_startr1, t_endr1, t_startr2, t_endr2, t_end, 
#         discrete, discrete_rate_timestep, tstops
#     )
# end