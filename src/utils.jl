"""
    tconvert(t, from_unit, to_unit)

Converts a time from one unit (`from_unit`) to another (`to_unit`).

Supported units (with accepted abbreviations) are:

* `picoseconds` (`ps`)
* `nanoseconds` (`ns`)
* `microseconds` (`us`)
* `milliseconds` (`ms`)
* `seconds` (`s`)
* `minutes` (`mins`)
* `hours` (`hrs`)
* `days`
* `months` (`mts`)
* `years` (`yrs`)

Will throw an error if unsupported units are provided.
"""
function tconvert(t::Real, from_unit::String, to_unit::String)
    if !(from_unit in keys(t_unit_map)) || !(to_unit in keys(t_unit_map))
        error("Unknown unit specified in time conversion!")
    end

    t = convert(Float64, t)
    t_conv = t * t_unit_map[from_unit] / t_unit_map[to_unit]

    return t_conv
end

"""
    tconvert(from_unit, to_unit)

Returns `tconvert(1.0, from_unit, to_unit)`.

Useful for just getting a conversion factor between time units, 
rather than converting a specific time directly.
"""
function tconvert(from_unit::String, to_unit::String)
    return tconvert(1.0, from_unit, to_unit)
end


"""
    tconvert(t, from_unit, to_unit)

Converts a vector of times from one unit (`from_unit`) to another (`to_unit`).

Supported units (with accepted abbreviations) are:

* `picoseconds` (`ps`)
* `nanoseconds` (`ns`)
* `microseconds` (`us`)
* `milliseconds` (`ms`)
* `seconds` (`s`)
* `minutes` (`mins`)
* `hours` (`hrs`)
* `days`
* `months` (`mts`)
* `years` (`yrs`)

Will throw an error if unsupported units are provided.
"""
function tconvert(t::Vector{<:Real}, from_unit::String, to_unit::String)
    if !(from_unit in keys(t_unit_map)) || !(to_unit in keys(t_unit_map))
        error("Unknown unit specified in time conversion!")
    end

    t = convert(Vector{Float64}, t)
    t_conv = t * t_unit_map[from_unit] / t_unit_map[to_unit]

    return t_conv
end


t_unit_map = Dict{String, Float64}(
    "picoseconds" => 1.0e-12,
    "ps" => 1.0e-12,
    "nanoseconds" => 1.0e-9,
    "ns" => 1.0e-9,
    "microseconds" => 1.0e-6,
    "us" => 1.0e-6,
    "milliseconds" => 1.0e-3,
    "ms" => 1.0e-3,
    "seconds" => 1.0,
    "s" => 1.0,
    "minutes" => 60.0,
    "mins" => 60.0,
    "hours" => 3600.0,
    "hrs" => 3600.0,
    "days" => 86400.0,
    "months" => 2.6297368e06,
    "mts" => 2.6297368e06,
    "years" => 3.15576e07,
    "yrs" => 3.15576e07
)


"""
    create_savepoints(start, stop, step)

Creates a range of savepoints, ensuring that the final time is included.

Safely modifies `step` to account for small floating point errors
introduced by calling `tconvert`.
"""
function create_savepoints(start::tType, stop::tType, step::tType) where {tType <: AbstractFloat}
    cstep = ((step > 1e-9) && (abs(step - floor(step)) < 1e-9)) ? round(step; sigdigits=9) : step
    r = collect(start:cstep:stop)
    if r[end] < stop
        push!(r, stop)
    end
    return r
end


"""
    trunc(Float64, x)

Overload of `Base.trunc` to allow for truncation of explicitly typed Real numbers.
"""
trunc(::Type{T}, x::T) where T<:Real = trunc(x)


"""
    DummyODEFunction{true}([syms])

Fake `ODEFunction` for use in `DummyODEProblem`.
"""
struct DummyODEFunction{iip} <: SciMLBase.AbstractSciMLFunction{iip} 
    syms
end

"""
    DummyODEProblem([uType, tType, u0, tspan, syms])

Fake `ODEProblem` implementing bare minimum fields for working with other SciMLBase code.

Needed within 'fake' `ODESolution`s for interpolation and plotting.

Symbolic names can be passed to the variables in the surrounding
`ODESolution` through the `syms` argument.
"""
struct DummyODEProblem{uType, tType, isinplace} <: SciMLBase.AbstractODEProblem{uType, tType, isinplace}
    u0
    tspan
    p
    f
end
function DummyODEProblem(; uType=Float64, tType=Float64, u0=[0.0], tspan=[0.0, 1.0], syms=nothing)
    return DummyODEProblem{uType, tType, true}(u0, tspan, nothing, DummyODEFunction{true}(syms))
end

