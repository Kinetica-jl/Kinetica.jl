"""
"""
struct ConditionSet{tType}
    symbols::Vector{Symbol}
    profiles::Vector{<:AbstractConditionProfile}
    discrete_updates::Bool
    ts_update::tType
end

"""
"""
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
    return ConditionSet(symbols, profiles, false, nothing)
end

"""
"""
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
    return ConditionSet(symbols, profiles, true, ts_update)
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

function get_profile(cs::ConditionSet, sym::Symbol)
    if !(sym in cs.symbols)
        throw(ErrorException("Condition $sym does not exist in this ConditionSet"))
    end
    loc = findfirst(==(sym), cs.symbols)
    return cs.profiles[loc]
end

function get_tstops(cs::ConditionSet)
    if isstatic(cs) throw(ErrorException("No tstops available, all conditions in ConditionSet are static.")) end
    all_tstops = [profile.tstops for profile in cs.profiles if isvariable(profile)]
    return sort(unique(reduce(vcat, all_tstops)))
end

function get_t_final(cs::ConditionSet)
    if isstatic(cs) throw(ErrorException("No t_end available, all conditions in ConditionSet are static.")) end
    all_t_end = [profile.t_end for profile in cs.profiles if isvariable(profile)]
    return maximum(all_t_end)
end