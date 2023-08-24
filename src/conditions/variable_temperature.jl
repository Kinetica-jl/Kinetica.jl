"""
Definitions for all variable temperature profile structs and theor outer constructor methods.

For compatibility, all temperature profile structs must implement the following fields:

* `Tgrad<:Function`
* `T_start<:AbstractFloat`
* `T_max<:AbstractFloat`
* `t_end<:AbstractFloat`
* `discrete_rates::Bool`
* `discrete_rate_timestep::Union{<:AbstractFloat, Nothing}`
* `tstops::Vector{<:AbstractFloat}`
"""


"""
Container for null temperature profile data and temperature gradient function.

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
* Timestep for discrete rate updates (`discrete_rate_timestep`)
* Times for the ODE solver to ensure calculation at (`tstops`)
"""
struct NullTprofile{uType, tType} <: AbstractTprofile
    Tgrad::Function
    T_start::uType
    T_max::uType
    t_end::tType
    discrete_rates::Bool
    discrete_rate_timestep::Union{tType, Nothing}
    tstops::Vector{tType}
end

"""
    Tprofile = NullTProfile(; T=T, t_end=t_end[, discrete_rate_timestep=discrete_rate_timestep])

Outer constructor for null temperature profile.

Should only be used for testing purposes (see struct
documentation).
"""
function NullTprofile(;
    T::uType,
    t_end::tType,
    discrete_rate_timestep::Union{tType, Nothing} = nothing
) where {uType <: AbstractFloat, tType <: AbstractFloat}
    function Tgrad(t)
        return 0.0
    end

    discrete = false
    if !isnothing(discrete_rate_timestep)
        if discrete_rate_timestep > t_end throw(ArgumentError("Error defining tstops, discrete_rate_timestep is too large.")) end
        discrete = true
        tstops = collect(0.0:discrete_rate_timestep:t_end)
    else
        tstops = [t_end]
    end

    return NullTprofile(Tgrad, T, T, t_end, discrete, discrete_rate_timestep, tstops)
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
* Timestep for discrete rate updates (`discrete_rate_timestep`)
* Times for the ODE solver to ensure calculation at (`tstops`)
"""
struct LinearTprofile{uType, tType} <: AbstractTprofile
    Tgrad::Function
    heatrate::uType
    T_start::uType
    T_end::uType
    T_max::uType
    t_end::tType
    discrete_rates::Bool
    discrete_rate_timestep::Union{tType, Nothing}
    tstops::Vector{tType}
end

"""
    Tprofile = LinearTprofile(; heatrate=heatrate, T_start=T_start, T_end=T_end
        discrete_rate_timestep=discrete_rate_timestep)

Outer constructor for linear temperature ramp temperature profile.

Determines the simulation end time from the provided temperatures
and gradient, then constructs the temperature gradient function 
(which returns `heatrate` for every timestep).

Also creates the `tstops` array for passing to adaptive timestepping
ODE solvers.
"""
function LinearTprofile(;
    heatrate::uType,
    T_start::uType,
    T_end::uType,
    discrete_rate_timestep::Union{tType, Nothing}=nothing
) where {uType <: AbstractFloat, tType <: AbstractFloat}
    t_end = (T_end - T_start)/heatrate

    function Tgrad(t)
        return heatrate
    end

    discrete = false
    if !isnothing(discrete_rate_timestep)
        t_end = tType(t_end)
        if discrete_rate_timestep > t_end throw(ArgumentError("Error defining tstops, discrete_rate_timestep is too large.")) end
        tstops = create_savepoints(tType(0.0), t_end, discrete_rate_timestep)
        discrete = true
    else
        tstops = [t_end]
    end

    return LinearTprofile(Tgrad, heatrate, T_start, T_end, T_end, t_end, discrete, discrete_rate_timestep, tstops)
end


"""
Container for simple double temperature ramp profile data and temperature gradient function.

This temperature profile steps through an initial plateau, 
followed by a discontinuous ramp to a middle temperature,
followed by another plateau, followed by another discontinuous
ramp to a final temperature, ending after a final plateau.

Note that while both `rate1` and `rate2` are defined as heating
rates, a negative value for either/both indicates a downward
temperature ramp, i.e. a cooling.

Contains fields for:
* Temperature gradient function (`Tgrad`)
* 1st heating rate (`rate1`)
* 2nd heating rate (`rate2`)
* Initial temperature (`T_start`)
* Middle plateau temperature (`T_mid`)
* Final temperature (`T_end`)
* Maximum temperature (`T_max`)
* Length of time to spend at initial temperature plateau (`t_startplat`)
* Length of time to spend at maximum temperature plateau (`t_midplat`)
* Length of time to spend at final temperature plateau (`t_endplat`)
* Time to start 1st heating ramp (`t_startr1`)
* Time to stop 1st heating ramp (`t_endr1`)
* Time to start 2nd heating ramp (`t_startr2`)
* Time to stop 2nd heating ramp (`t_endr2`)
* Time to stop calculation (`t_end`)
* Whether discrete rate updates should be enabled (`discrete_rates`)
* Timestep for discrete rate updates (`discrete_rate_timestep`)
* Times for the ODE solver to ensure calculation at (`tstops`)
"""
struct SimpleDoubleRampTprofile{uType, tType} <: AbstractTprofile
    Tgrad::Function
    rate1::uType
    rate2::uType
    T_start::uType
    T_mid::uType
    T_end::uType
    T_max::uType
    t_startplat::tType
    t_midplat::tType
    t_endplat::tType
    t_startr1::tType
    t_endr1::tType
    t_startr2::tType
    t_endr2::tType
    t_end::tType
    discrete_rates::Bool
    discrete_rate_timestep::Union{tType, Nothing}
    tstops::Vector{tType}
end

"""
    Tprofile = SimpleDoubleRampTprofile(; rate1=rate1, rate2=rate2, 
        T_start=T_start, T_mid=T_mid, T_end=T_end, t_startplat=t_startplat,
        t_midplat=t_midplat, t_endplat=t_endplat[, discrete_rate_timestep=discrete_rate_timestep])

Outer constructor for simple double temperature ramp temperature profile.

Creates all of the required absolute times from the provided 
relative times, then constructs the temperature gradient function
for the specified simple temperature ramp profile.

Also creates the `tstops` array for passing to adaptive timestepping
ODE solvers.
"""
function SimpleDoubleRampTprofile(;
    rate1::uType,
    rate2::uType,
    T_start::uType,
    T_mid::uType,
    T_end::uType,
    t_startplat::tType,
    t_midplat::tType,
    t_endplat::tType,
    discrete_rate_timestep::Union{tType, Nothing} = nothing
) where {uType <: AbstractFloat, tType <: AbstractFloat}
    if (T_mid < T_start && rate1 > 0) || (T_mid > T_end && rate2 < 0)
        error("Impossible temperature ramps defined. Check heating rates have the correct signs.")
    end

    t_startr1 = t_startplat
    t_endr1 = tType(t_startr1+ ((T_mid - T_start)/rate1))
    t_startr2 = t_endr1 + t_midplat
    t_endr2 = tType(t_startr2 + ((T_end - T_mid)/rate2))
    t_end = t_endr2 + t_endplat
    T_max = maximum([T_start, T_mid, T_end])

    function Tgrad(t)
        return uType(
            ((t < t_startr1) * 0.0) +
            ((t >= t_startr1 && t < t_endr1) * rate1) +
            ((t >= t_endr1 && t < t_startr2) * 0.0) +
            ((t >= t_startr2 && t < t_endr2) * rate2) +
            ((t >= t_endr2) * 0.0) 
        )
    end

    discrete = false
    if !isnothing(discrete_rate_timestep)
        if discrete_rate_timestep > t_end throw(ArgumentError("Error defining tstops, discrete_rate_timestep is too large.")) end
        tstops = reduce(vcat,
            [
                [tType(0.0)],
                create_savepoints(t_startr1, t_endr1, discrete_rate_timestep),
                create_savepoints(t_startr2, t_endr2, discrete_rate_timestep),
                [t_end]
            ])
        discrete = true
    else
        tstops = [
            t_startr1, t_endr1,
            t_startr2, t_endr2,
            t_end
        ]
    end

    return SimpleDoubleRampTprofile(
        Tgrad, rate1, rate2, T_start, T_mid, T_end, T_max, t_startplat, t_midplat,
        t_endplat, t_startr1, t_endr1, t_startr2, t_endr2, t_end, discrete,
        discrete_rate_timestep, tstops
    )
end


"""
Container for smooth double temperature ramp profile data and temperature gradient function.

This temperature profile steps through an initial plateau, 
followed by a smooth ramp up to a maximum temperature,
followed by another plateau, followed by a smooth ramp down
to a final temperature, ending after a final plateau.

Note that while both `rate1` and `rate2` are defined as heating
rates, a negative value for either/both indicates a downward
temperature ramp, i.e. a cooling.

Contains fields for:
* Temperature gradient function (`Tgrad`)
* 1st heating rate (`rate1`)
* 2nd heating rate (`rate2`)
* Initial temperature (`T_start`)
* Middle plateau temperature (`T_mid`)
* Final temperature (`T_end`)
* Maximum temperature (`T_max`)
* Length of time to spend at initial temperature plateau (`t_startplat`)
* Length of time to spend at maximum temperature plateau (`t_midplat`)
* Length of time to spend at final temperature plateau (`t_endplat`)
* Blending time to be taken out of each plateau (`t_blend`)
* Time to start 1st heating ramp (`t_startr1`)
* Time to stop 1st heating ramp (`t_endr1`)
* Time to start 2nd heating ramp (`t_startr2`)
* Time to stop 2nd heating ramp (`t_endr2`)
* Time to stop calculation (`t_end`)
* Whether discrete rate updates should be enabled (`discrete_rates`)
* Timestep for discrete rate updates (`discrete_rate_timestep`)
* Times for the ODE solver to ensure calculation at (`tstops`)
"""
struct SmoothDoubleRampTprofile{uType,tType} <: AbstractTprofile
    Tgrad::Function
    rate1::uType
    rate2::uType
    T_start::uType
    T_mid::uType
    T_end::uType
    T_max::uType
    t_startplat::tType
    t_midplat::tType
    t_endplat::tType
    t_blend::tType
    t_startr1::tType
    t_endr1::tType
    t_startr2::tType
    t_endr2::tType
    t_end::tType
    discrete_rates::Bool
    discrete_rate_timestep::Union{tType, Nothing}
    tstops::Vector{tType}
end

"""
    Tprofile = SmoothDoubleRampTprofile(; heatrate=heatrate, coolrate=coolrate, 
        T_start=T_start, T_max=T_max, T_end=T_end, t_startplateau=t_startplateau,
        t_midplateau=t_midplateau, t_endplateau=t_endplateau, t_blend=t_blend
        [, discrete_rate_timestep=discrete_rate_timestep])

Outer constructor for smooth double temperature ramp temperature profile.

Creates all of the required absolute times from the provided 
relative times, then constructs the temperature gradient function
for the specified smooth temperature ramp profile.

Also creates the `tstops` array for passing to adaptive timestepping
ODE solvers.
"""
function SmoothDoubleRampTprofile(;
    rate1::uType,
    rate2::uType,
    T_start::uType,
    T_mid::uType,
    T_end::uType,
    t_startplat::tType,
    t_midplat::tType,
    t_endplat::tType,
    t_blend::tType,
    discrete_rate_timestep::Union{tType, Nothing} = nothing
) where {uType <: AbstractFloat, tType <: AbstractFloat}
    if (T_mid < T_start && rate1 > 0) || (T_mid > T_end && rate2 < 0)
        error("Impossible temperature ramps defined. Check heating rates have the correct signs.")
    end
    if t_blend < 0.0
        error("t_blend must be a positive value.")
    end

    t_startr1 = t_startplat
    t_endr1 = tType(t_startr1+ ((T_mid - T_start)/rate1))
    t_startr2 = t_endr1 + t_midplat
    t_endr2 = tType(t_startr2 + ((T_end - T_mid)/rate2))
    t_end = t_endr2 + t_endplat
    T_max = maximum([T_start, T_mid, T_end])

    function Tgrad(t)
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

    discrete = false
    if !isnothing(discrete_rate_timestep)
        if discrete_rate_timestep > t_end throw(ArgumentError("Error defining tstops, discrete_rate_timestep is too large.")) end
        tstops = reduce(vcat,
            [
                [tType(0.0)],
                create_savepoints(t_startr1-t_blend, t_endr1+t_blend, discrete_rate_timestep),
                create_savepoints(t_startr2-t_blend, t_endr2+t_blend, discrete_rate_timestep),
                [t_end]
            ])
        discrete = true
    else
        tstops = [
            t_startr1-t_blend, t_startr1, t_startr1+t_blend, 
            t_endr1-t_blend, t_endr1, t_endr1+t_blend,
            t_startr2-t_blend, t_startr2, t_startr2+t_blend,
            t_endr2-t_blend, t_endr2, t_endr2+t_blend
        ]
    end

    return SmoothDoubleRampTprofile(
        Tgrad, rate1, rate2, T_start, T_mid, T_end, T_max, t_startplat, t_midplat,
        t_endplat, t_blend, t_startr1, t_endr1, t_startr2, t_endr2, t_end, 
        discrete, discrete_rate_timestep, tstops
    )
end