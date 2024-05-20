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
        [, ts_update=nothing]
    ))

Container for all conditions in a kinetic simulation.

Conditions can be static or variable, and variable
conditions can be gradient-based or directly usable.

Contains fields for:
* Symbolic representation of conditions (`symbols`)
* Condition profile for each symbol (`profiles`)
* Whether discrete rate constant updates are enabled for the conditions in this condition set (`discrete_updates`)
* Discrete rate constant update timestep, is `nothing` if `discrete_updates = false` (`ts_update`)

Constructor separates condition profiles from their symbols 
and parses numeric profiles into `StaticConditionProfile`s. 
Registers `AbstractVariableProfile`s with Symbolics to allow
for proper computation down the chain.

If `ts_update` is provided as a keyword argument, creates 
`tstops` arrays within each variable profile for use in 
discrete rate update simulations.
"""
function ConditionSet(d::Dict{Symbol, <:Any}; 
                      ts_update::Union{tType, Nothing}=nothing) where {tType <: AbstractFloat}
    symbols = collect(keys(d))
    profiles = AbstractConditionProfile[]
    for sym in symbols
        if d[sym] isa Number
            push!(profiles, StaticConditionProfile(d[sym]))
        elseif d[sym] isa AbstractConditionProfile
            if !isnothing(ts_update) create_discrete_tstops!(d[sym], ts_update) end
            push!(profiles, d[sym])
        else
            throw(ArgumentError("Condition $(sym) does not have a valid profile."))
        end
    end

    if isnothing(ts_update)
        cs = ConditionSet(symbols, profiles, false, nothing)
    else
        cs = ConditionSet(symbols, profiles, true, ts_update)
    end
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
    get_profile(cs::ConditionSet, sym::Symbol)

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
    get_initial_conditions(conditions::ConditionSet)

Extract initial values of conditions from `conditions`.

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
    get_static_conditions(conditions::ConditionSet)

Extracts static conditions from `conditions`.

Returns an array of `Pair`s linking Symbols to
static values. 
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
    get_variable_conditions(conditions::ConditionSet)

Extracts `ODESolution`s of variable conditions from `conditions`.

Returns an array of `Pair`s linking Symbols to
`ODESolution`s. 
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
    get_tstops(cs::ConditionSet)

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
    get_t_final(cs::ConditionSet)

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
    register_direct_conditions(cs::ConditionSet)

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
            pType = typeof(profile)
            eval(:(@register_symbolic($profile.f(t, p::$pType))))
        end
    end
end

"""
    register_gradient_conditions(cs::ConditionSet)

Calls `@register_symbolic` on all `AbstractGradientProfile` functions.

Workaround for otherwise needing to call the macro in the scope
of the `Main` module, once for each gradient profile. Instead,
handles this nicely without the user having to worry about explicitly
adding runtime-generated functions to the computation graph.
"""
function register_gradient_conditions(cs::ConditionSet)
    for sym in cs.symbols
        profile = get_profile(cs, sym)
        if isvariable(profile) && isgradientprofile(profile)
            pType = typeof(profile)
            eval(:(@register_symbolic($profile.grad(t, p::$pType))))
        end
    end
end


"""
    solve_variable_conditions!(cs::ConditionSet, pars::ODESimulationParams[, reset=false, solver=OwrenZen5(), solve_kwargs])

Solves all variable condition profiles over the timespan in `pars.tspan`.

Places all condition profile solutions in their `sol` field. In the
case of `AbstractDirectProfile`s, this creates a `DiffEqArray` to mimic
the regular DiffEq solver interface.

If condition profiles already exist, they will not be overwritten unless
`reset=true`. 

The `solver` and `solve_kwargs` arguments are used when solving gradient
profiles. The OwrenZen5 solver has shown to be a stable, accurate solver
capable of handling sudden gradient changes, and is a sensible default.
`solve_kwargs` is a `Dict` of keyword arguments that get passed to the 
`solve` call, with default state:

    solve_kwargs=Dict{Symbol, Any}(
        :abstol => 1e-6,
        :reltol => 1e-4
    )

These have been shown to be sensible defaults for most gradient profiles.
"""
function solve_variable_conditions!(cs::ConditionSet, pars::ODESimulationParams; 
    reset=false, solver=OwrenZen5(), solve_kwargs=Dict{Symbol, Any}(:abstol => 1e-6, :reltol => 1e-4))

    for (sym, profile) in zip(cs.symbols, cs.profiles)
        if isvariable(profile)
            solve_variable_condition!(profile, pars; sym, reset, solver, solve_kwargs)
        end
    end
end