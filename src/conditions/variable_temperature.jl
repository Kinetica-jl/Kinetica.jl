"""
Container for null temperature ramp profile data and temperature gradient function.

This temperature profile should only be used for debugging,
as it has a temperature gradient function which always
returns zero (i.e. there is no temperature change from
`T_start`). If constant temperature is required, the
regular `ODESolve` should always be used instead of an
`ODETemperatureSolve` with this temperature profile.

Contains fields for:
* Temperature gradient function (`Tgrad`)
* Initial temperature (`T_start`)
* Maximum plateau temperature, provided for internal consistency (`T_max`)
* Time to stop calculation (`t_end`)
* Whether discrete rate updates should be enabled (`discrete_rates`)
* Timestep for discrete rate updates (`rate_timestep`)
* Times for the ODE solver to ensure calculation at (`tstops`)
"""
struct NullTemperatureProfile{uType, tType} <: AbstractTemperatureProfile
    Tgrad::Function
    T_start::uType
    T_max::uType
    t_end::tType
    discrete_rates::Bool
    rate_timestep::Union{tType, Nothing}
    tstops::Vector{tType}
end

"""
    Tprofile = NullTemperatureProfile(; T=T, t_end=t_end[, rate_timestep=rate_timestep])

Outer constructor for null temperature profile.

Should only be used for testing purposes (see struct
documentation).
"""
function NullTemperatureProfile(;
    T::uType,
    t_end::tType,
    rate_timestep::Union{tType, Nothing} = nothing
) where {uType <: AbstractFloat, tType <: AbstractFloat}
    function Tgrad(t)
        return 0.0
    end

    discrete = false
    if !isnothing(rate_timestep)
        if rate_timestep > t_end throw(ArgumentError("Error defining tstops, rate_timestep is too large.")) end
        discrete = true
        tstops = collect(0.0:rate_timestep:t_end)
    else
        tstops = [t_end]
    end

    return NullTemperatureProfile(Tgrad, T, T, t_end, discrete, rate_timestep, tstops)
end


"""
Container for linear temperature ramp profile data and temperature gradient function.

This temperature profile represents a linear temperature
increase/decrease from `T_start` to `T_end`.

Contains fields for:
* Temperature gradient function (`Tgrad`)
* Heating rate, can be negative for cooling (`heatrate`)
* Initial temperature (`T_start`)
* Final temperature (`T_end`)
* Maximum temperature, used for rate cutoffs (`T_max`)
* Time to stop calculation (`t_end`)
* Whether discrete rate updates should be enabled (`discrete_rates`)
* Timestep for discrete rate updates (`rate_timestep`)
* Times for the ODE solver to ensure calculation at (`tstops`)
"""
struct LinearTemperatureProfile{uType, tType} <: AbstractTemperatureProfile
    Tgrad::Function
    heatrate::uType
    T_start::uType
    T_end::uType
    T_max::uType
    t_end::tType
    discrete_rates::Bool
    rate_timestep::Union{tType, Nothing}
    tstops::Vector{tType}
end

"""
    Tprofile = LinearTemperatureProfile(; heatrate=heatrate, T_start=T_start, T_end=T_end
        rate_timestep=rate_timestep)

Outer constructor for linear temperature ramp temperature profile.

Determines the simulation end time from the provided temperatures
and gradient, then constructs the temperature gradient function 
(which returns `heatrate` for every timestep).

Also creates the `tstops` array for passing to adaptive timestepping
ODE solvers.
"""
function LinearTemperatureProfile(;
    heatrate::uType,
    T_start::uType,
    T_end::uType,
    rate_timestep::Union{tType, Nothing}=nothing
) where {uType <: AbstractFloat, tType <: AbstractFloat}
    t_end = (T_end - T_start)/heatrate

    function Tgrad(t)
        return heatrate
    end

    discrete = false
    if !isnothing(rate_timestep)
        t_end = tType(t_end)
        if rate_timestep > t_end throw(ArgumentError("Error defining tstops, rate_timestep is too large.")) end
        tstops = create_savepoints(tType(0.0), t_end, rate_timestep)
        discrete = true
    else
        tstops = [t_end]
    end

    return LinearTemperatureProfile(Tgrad, heatrate, T_start, T_end, T_end, t_end, discrete, rate_timestep, tstops)
end


"""
Container for simple temperature ramp profile data and temperature gradient function.

This temperature profile steps through an initial plateau, 
followed by a discontinuous ramp up to a maximum temperature,
followed by another plateau, followed by another discontinuous
ramp down to a final temperature, ending after a final plateau.

Contains fields for:
* Temperature gradient function (`Tgrad`)
* Heating rate (`heatrate`)
* Cooling rate (`coolrate`)
* Initial temperature (`T_start`)
* Maximum plateau temperature (`T_max`)
* Final temperature (`T_end`)
* Length of time to spend at initial temperature plateau (`t_startplateau`)
* Length of time to spend at maximum temperature plateau (`t_midplateau`)
* Length of time to spend at final temperature plateau (`t_endplateau`)
* Time to start heating (`t_heatstart`)
* Time to stop heating (`t_heatend`)
* Time to start cooling (`t_coolstart`)
* Time to stop cooling (`t_coolend`)
* Time to stop calculation (`t_end`)
* Whether discrete rate updates should be enabled (`discrete_rates`)
* Timestep for discrete rate updates (`rate_timestep`)
* Times for the ODE solver to ensure calculation at (`tstops`)
"""
struct SimpleUpDownTemperatureProfile{uType, tType} <: AbstractTemperatureProfile
    Tgrad::Function
    heatrate::uType
    coolrate::uType
    T_start::uType
    T_max::uType
    T_end::uType
    t_startplateau::tType
    t_midplateau::tType
    t_endplateau::tType
    t_heatstart::tType
    t_heatend::tType
    t_coolstart::tType
    t_coolend::tType
    t_end::tType
    discrete_rates::Bool
    rate_timestep::Union{tType, Nothing}
    tstops::Vector{tType}
end

"""
    Tprofile = SimpleUpDownTemperatureProfile(; heatrate=heatrate, coolrate=coolrate, 
        T_start=T_start, T_max=T_max, T_end=T_end, t_startplateau=t_startplateau,
        t_midplateau=t_midplateau, t_endplateau=t_endplateau[, rate_timestep=rate_timestep])

Outer constructor for simple temperature ramp temperature profile.

Creates all of the required absolute times from the provided 
relative times, then constructs the temperature gradient function
for the specified simple temperature ramp profile.

Also creates the `tstops` array for passing to adaptive timestepping
ODE solvers.
"""
function SimpleUpDownTemperatureProfile(;
    heatrate::uType,
    coolrate::uType,
    T_start::uType,
    T_max::uType,
    T_end::uType,
    t_startplateau::tType,
    t_midplateau::tType,
    t_endplateau::tType,
    rate_timestep::Union{tType, Nothing} = nothing
) where {uType <: AbstractFloat, tType <: AbstractFloat}
    t_heatstart = t_startplateau
    t_heatend = tType(t_startplateau + ((T_max - T_start)/heatrate))
    t_coolstart = t_heatend + t_midplateau
    t_coolend = tType(t_coolstart + ((T_end - T_max)/coolrate))
    t_end = t_coolend + t_endplateau

    function Tgrad(t)
        return uType(
            ((t < t_heatstart) * 0.0) +
            ((t >= t_heatstart && t < t_heatend) * heatrate) +
            ((t >= t_heatend && t < t_coolstart) * 0.0) +
            ((t >= t_coolstart && t < t_coolend) * coolrate) +
            ((t >= t_coolend) * 0.0) 
        )
    end

    discrete = false
    if !isnothing(rate_timestep)
        if rate_timestep > t_end throw(ArgumentError("Error defining tstops, rate_timestep is too large.")) end
        tstops = reduce(vcat,
            [
                [tType(0.0)],
                create_savepoints(t_heatstart, t_heatend, rate_timestep),
                create_savepoints(t_coolstart, t_coolend, rate_timestep),
                [t_end]
            ])
        discrete = true
    else
        tstops = [
            t_heatstart, t_heatend,
            t_coolstart, t_coolend,
            t_end
        ]
    end

    return SimpleUpDownTemperatureProfile(
        Tgrad, heatrate, coolrate, T_start, T_max, T_end, t_startplateau, t_midplateau,
        t_endplateau, t_heatstart, t_heatend, t_coolstart, t_coolend, t_end, rate_timestep,
        discrete, tstops
    )
end


"""
Container for smooth temperature ramp profile data and temperature gradient function.

This temperature profile steps through an initial plateau, 
followed by a smooth ramp up to a maximum temperature,
followed by another plateau, followed by a smooth ramp down
to a final temperature, ending after a final plateau.

Contains fields for:
* Temperature gradient function (`Tgrad`)
* Heating rate (`heatrate`)
* Cooling rate (`coolrate`)
* Initial temperature (`T_start`)
* Maximum plateau temperature (`T_max`)
* Final temperature (`T_end`)
* Length of time to spend at initial temperature plateau (`t_startplateau`)
* Length of time to spend at maximum temperature plateau (`t_midplateau`)
* Length of time to spend at final temperature plateau (`t_endplateau`)
* Blending time to be taken out of each plateau (`t_blend`)
* Time to start heating (`t_heatstart`)
* Time to stop heating (`t_heatend`)
* Time to start cooling (`t_coolstart`)
* Time to stop cooling (`t_coolend`)
* Time to stop calculation (`t_end`)
* Whether discrete rate updates should be enabled (`discrete_rates`)
* Timestep for discrete rate updates (`rate_timestep`)
* Times for the ODE solver to ensure calculation at (`tstops`)
"""
struct SmoothUpDownTemperatureProfile{uType,tType} <: AbstractTemperatureProfile
    Tgrad::Function
    heatrate::uType
    coolrate::uType
    T_start::uType
    T_max::uType
    T_end::uType
    t_startplateau::tType
    t_midplateau::tType
    t_endplateau::tType
    t_blend::tType
    t_heatstart::tType
    t_heatend::tType
    t_coolstart::tType
    t_coolend::tType
    t_end::tType
    discrete_rates::Bool
    rate_timestep::Union{tType, Nothing}
    tstops::Vector{tType}
end

"""
    Tprofile = SmoothUpDownTemperatureProfile(; heatrate=heatrate, coolrate=coolrate, 
        T_start=T_start, T_max=T_max, T_end=T_end, t_startplateau=t_startplateau,
        t_midplateau=t_midplateau, t_endplateau=t_endplateau, t_blend=t_blend
        [, rate_timestep=rate_timestep])

Outer constructor for smooth temperature ramp temperature profile.

Creates all of the required absolute times from the provided 
relative times, then constructs the temperature gradient function
for the specified smooth temperature ramp profile.

Also creates the `tstops` array for passing to adaptive timestepping
ODE solvers.
"""
function SmoothUpDownTemperatureProfile(;
    heatrate::uType,
    coolrate::uType,
    T_start::uType,
    T_max::uType,
    T_end::uType,
    t_startplateau::tType,
    t_midplateau::tType,
    t_endplateau::tType,
    t_blend::tType,
    rate_timestep::Union{tType, Nothing} = nothing
) where {uType <: AbstractFloat, tType <: AbstractFloat}
    t_heatstart = t_startplateau
    t_heatend = tType(t_startplateau + ((T_max - T_start)/heatrate))
    t_coolstart = t_heatend + t_midplateau
    t_coolend = tType(t_coolstart + ((T_end - T_max)/coolrate))
    t_end = t_coolend + t_endplateau

    function Tgrad(t)
        return uType(
            ((t < t_heatstart-t_blend) * 0.0) +
            ((t >= t_heatstart-t_blend && t < t_heatstart+t_blend) * (heatrate*(t-t_heatstart-t_blend)/(2*t_blend) + heatrate)) +
            ((t >= t_heatstart+t_blend && t < t_heatend-t_blend) * heatrate) +
            ((t >= t_heatend-t_blend && t < t_heatend+t_blend) * (-heatrate*(t-t_heatend-t_blend)/(2*t_blend))) +
            ((t >= t_heatend+t_blend && t < t_coolstart-t_blend) * 0.0) +
            ((t >= t_coolstart-t_blend && t < t_coolstart+t_blend) * (coolrate*(t-t_coolstart-t_blend)/(2*t_blend) + coolrate)) +
            ((t >= t_coolstart+t_blend && t < t_coolend-t_blend) * coolrate) +
            ((t >= t_coolend-t_blend && t < t_coolend+t_blend) * (-coolrate*(t-t_coolend-t_blend)/(2*t_blend))) +
            ((t >= t_coolend+t_blend) * 0.0)
        )
    end

    discrete = false
    if !isnothing(rate_timestep)
        if rate_timestep > t_end throw(ArgumentError("Error defining tstops, rate_timestep is too large.")) end
        tstops = reduce(vcat,
            [
                [tType(0.0)],
                create_savepoints(t_heatstart-t_blend, t_heatend+t_blend, rate_timestep),
                create_savepoints(t_coolstart-t_blend, t_coolend+t_blend, rate_timestep),
                [t_end]
            ])
        discrete = true
    else
        tstops = [
            t_heatstart-t_blend, t_heatstart, t_heatstart+t_blend, 
            t_heatend-t_blend, t_heatend, t_heatend+t_blend,
            t_coolstart-t_blend, t_coolstart, t_coolstart+t_blend,
            t_coolend-t_blend, t_coolend, t_coolend+t_blend
        ]
    end

    return SmoothUpDownTemperatureProfile(
        Tgrad, heatrate, coolrate, T_start, T_max, T_end, t_startplateau, t_midplateau,
        t_endplateau, t_blend, t_heatstart, t_heatend, t_coolstart, t_coolend, t_end, 
        rate_timestep, discrete, tstops
    )
end