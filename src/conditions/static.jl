"""
Definition for static condition profile.

Conditions defined this way are static for the duration
of a simulation.
"""
struct StaticConditionProfile{uType} <: AbstractStaticProfile
    value::uType
end

