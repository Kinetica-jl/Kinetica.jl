"""
Container for all conditions in a kinetic simulation.

Conditions can be static or variable, and variable
conditions can be gradient-based or directly usable.

Contains fields for:
* Symbolic representation of conditions (`symbols`)
* Condition profile for each symbol (`profiles`)
* Whether discrete rate constant updates are enabled for the conditions in this condition set (`discrete_updates`)
* Discrete rate constant update timestep, is `nothing` if `discrete_updates = false` (`ts_update`)
"""
struct ConditionSet{tType}
    symbols::Vector{Symbol}
    profiles::Vector{<:AbstractConditionProfile}
    discrete_updates::Bool
    ts_update::tType
end

"""
    conditions = ConditionSet(Dict(
        :C1 => ConditionType1(...),
        :C2 => ConditionType2(...))
        [, ts_update]
    ))

Outer constructor for `ConditionSet`..

Separates condition profiles from their symbols and parses
numeric profiles into `StaticConditionProfile`s. Registers
`AbstractDirectProfile`s with Symbolics to allow for proper
computation down the chain.

If `ts_update` is provided, creates `tstops` arrays within
each variable profile for use in discrete rate update 
simulations.
"""
function ConditionSet end

function ConditionSet(d::Dict{Symbol, <:Any})
    symbols = collect(keys(d))
    profiles = AbstractConditionProfile[]
    for sym in symbols
        if d[sym] isa Number
            push!(profiles, StaticConditionProfile(d[sym]))
        elseif d[sym] isa AbstractConditionProfile
            push!(profiles, d[sym])
        else
            throw(ArgumentError("Condition $(symbols[i]) does not have a valid profile."))
        end
    end
    cs = ConditionSet(symbols, profiles, false, nothing)
    register_direct_conditions(cs)
    register_gradient_conditions(cs)
    return cs
end

function ConditionSet(d::Dict{Symbol, <:Any}, ts_update::AbstractFloat)
    symbols = collect(keys(d))
    profiles = AbstractConditionProfile[]
    for sym in symbols
        if d[sym] isa Number
            push!(profiles, StaticConditionProfile(d[sym]))
        elseif d[sym] isa AbstractConditionProfile
            create_discrete_tstops(d[sym], ts_update)
            push!(profiles, d[sym])
        else
            throw(ArgumentError("Condition $(symbols[i]) does not have a valid profile."))
        end
    end
    cs = ConditionSet(symbols, profiles, true, ts_update)
    register_direct_conditions(cs)
    register_gradient_conditions(cs)
    return cs
end


function isstatic(cs::ConditionSet, sym::Symbol)
    if !(sym in cs.symbols)
        throw(ErrorException("Condition $sym does not exist in this ConditionSet"))
    end
    loc = findfirst(==(sym), cs.symbols)
    return isstatic(cs.profiles[loc])
end

function isstatic(cs::ConditionSet)
    return all([isstatic(profile) for profile in cs.profiles])
end

function isvariable(cs::ConditionSet, sym::Symbol)
    if !(sym in cs.symbols)
        throw(ErrorException("Condition $sym does not exist in this ConditionSet"))
    end
    loc = findfirst(==(sym), cs.symbols)
    return isvariable(cs.profiles[loc])
end

function isvariable(cs::ConditionSet)
    return all([isvariable(profile) for profile in cs.profiles])
end


"""
    profile = get_profile(cs, sym)

Gets the condition profile linked to `sym` from `cs`.
"""
function get_profile(cs::ConditionSet, sym::Symbol)
    if !(sym in cs.symbols)
        throw(ErrorException("Condition $sym does not exist in this ConditionSet"))
    end
    loc = findfirst(==(sym), cs.symbols)
    return cs.profiles[loc]
end


"""
    initial_conditions = get_initial_conditions(conditions)

Extract initial conditions from ConditionSet.

Returns an array of `Pair`s linking Symbols to
initial values. For `AbstractStaticProfile`s,
initial values are their static values. For
`AbstractVariableProfile`s, initial values are
their `X_start` values.
"""
function get_initial_conditions(conditions::ConditionSet)
    ics = []
    for (sym, prof) in zip(conditions.symbols, conditions.profiles)
        if isstatic(prof)
            push!(ics, Pair(sym, prof.value))
        else
            push!(ics, Pair(sym, prof.X_start))
        end
    end
    return ics
end


"""
    static_conditions = get_static_conditions(conditions)

Extract static conditions from ConditionSet.
"""
function get_static_conditions(conditions::ConditionSet)
    scs = []
    for (sym, prof) in zip(conditions.symbols, conditions.profiles)
        if isstatic(prof)
            push!(scs, Pair(sym, prof.value))
        end
    end
    return scs
end


"""
    variable_conditions = get_variable_conditions(conditions)

Extracts `ODESolution`s of variable conditions from ConditionSet.
"""
function get_variable_conditions(conditions::ConditionSet)
    vcs = []
    for (sym, prof) in zip(conditions.symbols, conditions.profiles)
        if isvariable(prof)
            push!(vcs, Pair(sym, prof.sol))
        end
    end
    return vcs
end


"""
    tstops = get_tstops(cs)

Retrieves a sorted array of unique time stops from all condition profiles in `cs`.

Should be used for passing a unified set of time stops
to a discrete rate constant update solver. Will throw an
error if all condition profiles are static, as they have
no `tstops`.
"""
function get_tstops(cs::ConditionSet)
    if isstatic(cs) throw(ErrorException("No tstops available, all conditions in ConditionSet are static.")) end
    all_tstops = [profile.tstops for profile in cs.profiles if isvariable(profile)]
    return sort(unique(reduce(vcat, all_tstops)))
end


"""
    t_final = get_t_final(cs)

Retrieves the last necessary time point needed to encompass all variable condition profiles in `cs`.

Will throw an error if all condition profiles are static,
they have no set endpoints.
"""
function get_t_final(cs::ConditionSet)
    if isstatic(cs) throw(ErrorException("No t_end available, all conditions in ConditionSet are static.")) end
    all_t_end = [profile.t_end for profile in cs.profiles if isvariable(profile)]
    return maximum(all_t_end)
end


"""
    register_direct_conditions(cs)

Calls `@register_symbolic` on all `AbstractDirectProfile` functions.

Workaround for otherwise needing to call the macro in the scope
of the `Main` module, once for each direct profile. Instead,
handles this nicely without the user having to worry about explicitly
adding runtime-generated functions to the computation graph.
"""
function register_direct_conditions(cs::ConditionSet)
    for sym in cs.symbols
        profile = get_profile(cs, sym)
        if isvariable(profile) && isdirectprofile(profile)
            eval(:(@register_symbolic($profile.f(t))))
        end
    end
end

"""
    register_gradient_conditions(cs)

Calls `@register_symbolic` on all `AbstractGradientProfile` functions.

Workaround for otherwise needing to call the macro in the scope
of the `Main` module, once for each direct profile. Instead,
handles this nicely without the user having to worry about explicitly
adding runtime-generated functions to the computation graph.
"""
function register_gradient_conditions(cs::ConditionSet)
    for sym in cs.symbols
        profile = get_profile(cs, sym)
        if isvariable(profile) && isgradientprofile(profile)
            eval(:(@register_symbolic($profile.grad(t))))
        end
    end
end


"""
    solve_variable_conditions!(cs, pars[, reset, kwargs...])

Solves all variable condition profiles over the timespan in `pars.tspan`.

Places all condition profile solutions in their `sol` field. In the
case of `AbstractDirectProfile`s, this creates an `ODESolution` to mimic
the regular DiffEq solver interface.
"""
function solve_variable_conditions!(cs::ConditionSet, pars::ODESimulationParams; reset=false, kwargs...)
    for profile in cs.profiles
        if isvariable(profile)
            solve_variable_condition!(profile, pars; reset, kwargs...)
        end
    end
end