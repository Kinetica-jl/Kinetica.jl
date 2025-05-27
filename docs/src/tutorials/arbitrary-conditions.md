# Arbitrary Simulation Conditions

```@setup arbitrary_conditions
using Kinetica, Plots
mkpath("../assets/tutorials/arbitrary_conditions")
cspars = ODESimulationParams(
    tspan = (0.0, 10.0),
    u0 = [0.0],
    solver = nothing
)
```

Alongside the modular kinetic calculator interface, the ability to pass arbitrary combinations of simulation conditions through to be used in CRN integration is one of the key elements of extensibility and customisation within Kinetica.

By only requiring that conditions use consistent symbolic names at their definition and within a kinetic calculator, any variable condition profile, defined either directly or by its gradient with respect to time, can be symbolically bound to any quantity of interest.

## [`ConditionSet`](@ref)

At the core of this system is the [`ConditionSet`](@ref), which acts as an aggregator for individual condition profiles and their symbols. At their definition, [`ConditionSet`](@ref)s take a dictionary of `Symbol => Profile()` mappings. Each `Symbol` can realistically be anything that Julia allows, but by convention (and for compatibility with most calculator implementations) we stick to the usual abbreviations for common conditions - `:T` for temperature, `:P` for pressure, `:V` for volume, etc. Each condition profile can be one of three options:

* A `Number` representing a static value for the given condition to take for the duration of the simulation. Internally this is converted into a [`Kinetica.StaticConditionProfile`](@ref), but this is just a container for the number within.
* A directly variable condition profile, e.g. [`LinearDirectProfile`](@ref). These are variable condition profiles where the condition is implemented directly as a function of time.
* A gradient-variable condition profile, e.g. [`LinearGradientProfile`](@ref). These are variable condition profiles where the condition is implemented indirectly through its gradient with respect to time. These profiles must be numerically integrated (handled automatically within Kinetica) before they can be used.

The [`ConditionSet`](@ref) constructor takes an optional keyword argument, `ts_update`. If provided, this argument causes any kinetic simulations done with this [`ConditionSet`](@ref) to use the discrete rate update approximation, which is usually desired in any moderate to large-scale CRN simulations. For more information, see the tutorial on [ODE Solution](@ref).

A [`ConditionSet`](@ref) which implements static simulation volume, linearly increasing temperature and linearly decreasing pressure with a rate constant update timestep of 1 ms could therefore look like the following:

```@example arbitrary_conditions
conditions = ConditionSet(Dict(
    :V => 1000.0,
    :T => LinearDirectProfile(;
        rate = 20.0,
        X_start = 300.0,
        X_end = 500.0
    ),
    :P => LinearGradientProfile(;
        rate = -50.0,
        X_start = 1e5,
        X_end = 9e4
    )),
    ts_update = tconvert(1.0, "ms", "s")
)
```

!!! note "Converting times"
    By default, Kinetica works in units of seconds (this can be changed within the kinetic calculator being used). While the conversion above may be a bit redundant, the [`tconvert`](@ref) function can be used to quickly convert between commonly used time units.

!!! warning "Discrete Rate Approximation Behaviour"
    The behaviour of the `ts_update` argument is subject to change in the near future. This will likely not be a dramatic change, but it may be worth bearing in mind.

### Useful Functions

Once constructed, [`ConditionSet`](@ref)s can be queried in a number of ways. To fetch any of the profiles within, [`get_profile`](@ref) can be called:

```@example arbitrary_conditions
get_profile(conditions, :T)
```

To test if a given condition is static or variable, the [`isstatic`](@ref) and [`isvariable`](@ref) functions can be called:

```@example arbitrary_conditions
println("Temperature profile is static: $(isstatic(conditions, :T))")
println("Pressure profile is variable: $(isvariable(conditions, :P))")
```

To get the final time at which all condition profiles have stopped varying, [`get_t_final`](@ref) can be called. This returns the maximum value of each condition profile's `t_end` attribute (see below):

```@example arbitrary_conditions
get_t_final(conditions)
```

## Condition Profile Showcase

Below are the currently implemented variable condition profiles, along with examples of their shapes. Condition profiles are being added as we need them, so this library is currently quite small. You can help us out by adding new profiles and submitting a pull request (see the Development section on [Condition Profiles](@ref)), or by requesting them to be added on our [Issues page](https://github.com/Kinetica-jl/Kinetica.jl/issues)!

### Directly Variable Condition Profiles

#### [`LinearDirectProfile`](@ref)

This profile represents a linear change from one value to another. It has a piecewise linear condition function, defined as follows for the arbitrary condition ``X``:

```math
X\left( t \right) = 
    \begin{cases}
        \texttt{X\_start}, & \text{if } t \leq 0.0 \\
        \texttt{X\_start} + t\left( \texttt{rate} \right), & \text{if } t > 0.0 \text{ and } t \leq t_{\text{end}} \\
        \texttt{X\_end}, & \text{if } t > t_{\text{end}} \\
    \end{cases}
```

where ``t_{\text{end}} = \left( \texttt{X\_end} - \texttt{X\_start} \right) / \texttt{rate}``. For example:

```@setup arbitrary_conditions
# This is not being run in an example block with lines hidden because there
# is seemingly no way of not having any kind of return block, and redefining
# the functions within profiles prints to stderr, which is returned when the
# example returns nothing.
cs = ConditionSet(Dict(
    :T => LinearDirectProfile(;
        X_start = 300.0,
        X_end = 500.0,
        rate=20.0)
))
cspars.tspan = (0.0, get_t_final(cs))
solve_variable_conditions!(cs, cspars)
plot(get_profile(cs, :T).sol, label="T", xlabel="t")
savefig("../assets/tutorials/arbitrary_conditions/lineardirectprofile.svg")
```

```julia
cs = ConditionSet(Dict(
    :T => LinearDirectProfile(;
        X_start = 300.0,
        X_end = 500.0,
        rate=20.0)
))
```

![](../assets/tutorials/arbitrary_conditions/lineardirectprofile.svg)

### Gradient-Variable Condition Profiles

#### [`LinearGradientProfile`](@ref)

The gradient-based implementation of [`LinearDirectProfile`](@ref). Either can be used, they should be equally accurate. Mostly serves as an example of how gradient profiles differ in implementation to their directly variable counterparts. It has a piecewise linear gradient function, defined as follows for the arbitrary condition ``X``:

```math
\frac{\mathrm{d} X\left( t \right)}{\mathrm{d}t} = 
    \begin{cases}
        0.0, & \text{if } t \leq 0.0 \\
        \texttt{rate}, & \text{if } t > 0.0 \text{ and } t \leq t_{\text{end}} \\
        0.0, & \text{if } t > t_{\text{end}} \\
    \end{cases}
```

where ``t_{\text{end}} = \left( \texttt{X\_end} - \texttt{X\_start} \right) / \texttt{rate}``. For example:

```@setup arbitrary_conditions
cs = ConditionSet(Dict(
    :P => LinearGradientProfile(;
        X_start = 1e5,
        X_end = 9e4,
        rate=-50.0)
))
cspars.tspan = (0.0, get_t_final(cs))
solve_variable_conditions!(cs, cspars)
plot(get_profile(cs, :P).sol)
savefig("../assets/tutorials/arbitrary_conditions/lineargradientprofile.svg")
```

```julia
cs = ConditionSet(Dict(
    :P => LinearGradientProfile(;
        X_start = 1e5,
        X_end = 9e4,
        rate=-50.0)
))
```

![](../assets/tutorials/arbitrary_conditions/lineargradientprofile.svg)

#### [`DoubleRampGradientProfile`](@ref)

This profile represents two linear condition ramps, each of which can have either a positive or a negative gradient, separated by a plateau of variable time. The profile also begins and ends with variable-length condition plateaus to enable equilibration at the initial and final values. It has a piecewise linear gradient function, defined as follows for the arbitrary condition ``X``:

```math
\frac{\mathrm{d} X\left( t \right)}{\mathrm{d}t} = 
    \begin{cases}
        0.0, & \text{if } t < t_{r1, \text{start}} \\
        \texttt{rate1}, & \text{if } t_{r1, \text{start}} \leq t < t_{r1, \text{end}} \\
        0.0, & \text{if } t_{r1, \text{end}} \leq t < t_{r2, \text{start}} \\
        \texttt{rate2}, & \text{if } t_{r2, \text{start}} \leq t < t_{r2, \text{end}} \\
        0.0, & \text{if } t \geq t_{r2, \text{end}} 
    \end{cases}
```

where ``\texttt{rate1}`` and ``\texttt{rate2}`` are the rates of change of the two linear ramps, and ``t_{r1, \text{start}}``, ``t_{r1, \text{end}}``, ``t_{r2, \text{start}}`` and ``t_{r2, \text{end}}`` are the respective start- and end-times of the first and second ramps, determined by the lengths of the starting, middle and ending plateaus.

The profile features an optional argument `t_blend`, which can be used to create smooth transitions between the otherwise discontinuous gradient changes through linear interpolation. For example:

```@example arbitrary_conditions
cs = ConditionSet(Dict(
    :K => DoubleRampGradientProfile(;
        X_start = 100.0,
        t_start_plateau = 3.0,
        rate1 = 10.0,
        X_mid = 250.0,
        t_mid_plateau = 5.0,
        rate2 = -25.0,
        X_end = 50.0,
        t_end_plateau = 10.0,
        t_blend = 0.5)
))
cspars.tspan = (0.0, get_t_final(cs)) # hide
solve_variable_conditions!(cs, cspars) # hide
plot(get_profile(cs, :K).sol) # hide
savefig("../assets/tutorials/arbitrary_conditions/doublerampgradientprofile.svg"); nothing # hide
```

![](../assets/tutorials/arbitrary_conditions/doublerampgradientprofile.svg)