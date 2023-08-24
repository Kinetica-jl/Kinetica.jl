abstract type AbstractConditionProfile end

# Single variable condition profiles
abstract type AbstractTprofile <: AbstractConditionProfile end
abstract type AbstractPprofile <: AbstractConditionProfile end
abstract type AbstractVprofile <: AbstractConditionProfile end

# Double variable condition profiles
abstract type AbstractTPprofile <: AbstractConditionProfile end
abstract type AbstractTVprofile <: AbstractConditionProfile end
abstract type AbstractPVprofile <: AbstractConditionProfile end

# Triple variable condition profiles
abstract type AbstractTPVprofile <: AbstractConditionProfile end