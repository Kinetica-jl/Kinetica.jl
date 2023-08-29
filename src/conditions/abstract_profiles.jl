abstract type AbstractConditionProfile end

abstract type AbstractStaticProfile <: AbstractConditionProfile end
abstract type AbstractVariableProfile <: AbstractConditionProfile end

abstract type AbstractGradientProfile <: AbstractVariableProfile end
abstract type AbstractDirectProfile <: AbstractVariableProfile end


function isstatic(profile::pType) where {pType <: AbstractStaticProfile}
    return true
end

function isstatic(profile::pType) where {pType <: AbstractVariableProfile}
    return false
end

function isvariable(profile::pType) where {pType <: AbstractStaticProfile}
    return false
end

function isvariable(profile::pType) where {pType <: AbstractVariableProfile}
    return true
end