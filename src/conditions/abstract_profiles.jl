abstract type AbstractConditionProfile end

abstract type AbstractStaticProfile <: AbstractConditionProfile end
abstract type AbstractVariableProfile <: AbstractConditionProfile end

abstract type AbstractGradientProfile <: AbstractVariableProfile end
abstract type AbstractDirectProfile <: AbstractVariableProfile end

"""
    isstatic(cs[, sym])

Determines if condition profiles in a `ConditionSet` are static.

When a Symbol `sym` is provided, only checks if the profile
linked to this Symbol is static. If no Symbol is provided,
checks is all profiles are static.

    isstatic(profile)

Determines if a given condition profile is static.
"""
function isstatic end

function isstatic(::pType) where {pType <: AbstractStaticProfile}
    return true
end

function isstatic(::pType) where {pType <: AbstractVariableProfile}
    return false
end

"""
    isvariable(cs[, sym])

Determines if condition profiles in a `ConditionSet` are variable.

When a Symbol `sym` is provided, only checks if the profile
linked to this Symbol is variable. If no Symbol is provided,
checks is all profiles are variable.

    isvariable(profile)

Determines if a given condition profile is variable.
"""
function isvariable end

function isvariable(::pType) where {pType <: AbstractStaticProfile}
    return false
end

function isvariable(::pType) where {pType <: AbstractVariableProfile}
    return true
end


"""
    isgradientprofile(profile)

Detemines if an `AbstractVariableProfile` is gradient-based.
"""
function isgradientprofile end

function isgradientprofile(::pType) where {pType <: AbstractGradientProfile}
    return true
end

function isgradientprofile(::pType) where {pType <: AbstractDirectProfile}
    return false
end


"""
    isdirectprofile(profile)

Determines if an `AbstractVariableProfile` has a direct equation.
"""
function isdirectprofile end

function isdirectprofile(::pType) where {pType <: AbstractGradientProfile}
    return false
end

function isdirectprofile(::pType) where {pType <: AbstractDirectProfile}
    return true    
end


"""
    create_discrete_tstops(profile, ts_update)

Creates a custom array of time stops within `profile.tstops`.

This array contains a time stop every `ts_update`, but attempts
to intellegently avoid unnecessary time stops in areas where they
are not needed, i.e. when the given profile is stationary.
"""
function create_discrete_tstops end


"""
    minimum(profile)

Determines the minimum value of an `AbstractVariableProfile` from its solution.

Requires that the profile's solution has been generated
and is available from `profile.sol`. 

Note that this will only return the minimum value of a 
condition over the timespan requested during the solution,
so this may not be the true minimum value of the condition
if the solution ends before it is reached.
"""
function Base.minimum(profile::AbstractVariableProfile)
    if isnothing(profile.sol)
        throw(ErrorException("Condition profile is missing a solution."))
    end
    return minimum(profile.sol)
end


"""
    maximum(profile)

Determines the maximum value of an `AbstractVariableProfile` from its solution.

Requires that the profile's solution has been generated
and is available from `profile.sol`. 

Note that this will only return the maximum value of a 
condition over the timespan requested during the solution,
so this may not be the true maximum value of the condition
if the solution ends before it is reached.
"""
function Base.maximum(profile::AbstractVariableProfile)
    if isnothing(profile.sol)
        throw(ErrorException("Condition profile is missing a solution."))
    end
    return maximum(profile.sol)
end