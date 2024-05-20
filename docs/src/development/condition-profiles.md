# Condition Profiles

```@setup condition_profiles
using Kinetica
prev_res = load_output("../my_CRN_out/direct_network_final.bson")
sd, rd = prev_res.sd, prev_res.rd

using Sundials
pars = ODESimulationParams(
    tspan = (0.0, 1.0),
    u0 = Dict("C" => 1.0),
    solver = CVODE_BDF(; linear_solver=:KLU),
    abstol = 1e-11,
    reltol = 1e-9
)

using BSON
calc_pars = BSON.load("../../../examples/getting_started/arrhenius_params.bson")
calc = PrecalculatedArrheniusCalculator(calc_pars[:Ea], calc_pars[:A]; k_max=1e12)

mkpath("../assets/tutorials/condition_profiles")
```

For an introduction into Kinetica's arbitrary simulation condition framework, see the tutorial on [Arbitrary Simulation Conditions](@ref). This page aims to explain the underlying implementation of Kinetica's variable condition profiles, such that users can easily extend Kinetica by adding their own profiles to represent new macroscopic state variable changes.

As noted previously, variable condition profiles can take on one of two types of definition which effect how they are represented and interacted with inside Kinetica's solvers:

* Directly variable condition profiles (subtypes of `Kinetica.AbstractDirectProfile`) encode a condition ``X`` as a continuous function of time ``X(t)``. The result of this function is then computed directly whenever it is needed at a given time ``t``.
* Gradient-variable condition profiles (subtypes of `Kinetica.AbstractGradientProfile`) encode a condition ``X`` within a gradient with respect to time ``\frac{dX}{dt}``. These condition profiles must be integrated with time using DifferentialEquations.jl to yield an interpolable `ODESolution` which can be queried to produce values of the condition at time ``t``.

Both types of variable condition profile are simple to implement within Kinetica, provided a few criteria are fulfilled:

* Both are implemented as structs with supertypes as above that indicate how Kinetica should handle them. These structs must have a number of fields, but are otherwise entirely customisable for whatever data needs to be held. These usually come with outer constructors that create variables that are needed within the profile's condition function, as well as binding a condition function to one of the struct's fields.
* A condition function (usually called `_f_[struct name]` for direct condition profiles or `_grad_[struct_name]` for gradient condition profiles) must exist such that it can be bound to a condition profile struct. This function takes two arguments - the time `t` and the `profile` that it will ultimately be bound to.
* A method of [`Kinetica.create_discrete_tstops!`](@ref) must be created for the new profile struct which tells Kinetica where rate constants should be updated in discrete rate constant update simulations.

We'll demonstrate this by implementing both directly variable and gradient-variable versions of a simple sinusoidal condition profile.

!!! note "Why not a functor?"
    The implementation of condition functions described above may look a bit strange. Variable condition profiles require condition functions that always take the same arguments for internal consistency, but are also capable of taking a potentially large number of runtime-defined parameters that are accessible from some kind of `Kinetica.AbstractVariableProfile`. Why then do we pass a condition function as a field of the profile struct during construction and then pass the struct *back* to that function instead of using a [functor](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects), which would seemingly fill both of these roles?

    While the functor approach would certainly be more convenient, it's currently incompatible with [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) (more specifically its function registration API), which Kinetica uses to programatically embed variable condition profiles into systems of ODEs. Instead, we are forced to use a regular function, which then takes its parameters from a condition profile struct such that its call signature is always `conditionfunc(t, profile::AbstractVariableProfile)`. This function is bound to its respective condition profile struct so that irrespective of its true name, it is always accessible through `profile.f` for direct profiles, and `profile.grad` for gradient profiles.

    Similarly, we could define condition functions within the outer constructors, allowing them to implicitly make use of profile variables defined within the constructor and eliminating the need to pass an instance of the profile struct as an argument. However, due to how functions are registered within Symbolics, this can lead to condition functions overwriting one another. The separate condition function approach is admittedly a bit cumbersome, but it is currently the best viable approach that allows for hands-off creation of [`ConditionSet`](@ref)s for users.

## Directly Variable Profile Implementation

When a condition profile can be represented exactly as a simple function of time, it is usually easiest to implement as a directly variable condition profile. In the case of our sinusoidal profile, this function of time will be a sinusoid:
```math
X(t)=A\sin\left( 2\pi f t \right)
```
where ``A`` is the amplitude and ``f`` is the number of oscillations per second.

We will begin by creating the struct for our condition profile, which we'll call `SinusoidDirectProfile`. Aside from the parameters above, we'll also need to implement fields for some basics which Kinetica expects in every directly variable profile:

* `f`: The function which will be bound in the struct's outer constructor.
* `X_start`: The initial value of the condition profile.
* `t_end`: The time at which the condition profile should stop varying. This can either be provided by the user or calculated within the outer constructor.
* `tstops`: Time points at which ODE solvers should ensure to stop at and recalculate. This array is useful for handling discontinuities in condition profiles, and is usually modified by [`Kinetica.create_discrete_tstops!`](@ref).
* `sol`: Profile solution over the requested timespan, usually created internally by calling [`Kinetica.solve_variable_condition!(::Kinetica.AbstractDirectProfile, ::ODESimulationParams)`](@ref).

For this profile, we will also implement a field `X_end` to signify the value of the condition once it stops varying. In some profiles, this is more easily calculated by the profile itself, e.g. as a function of time when `t_end` is given. Conversely, in other profiles, it can be convenient to only provide `X_end` and let `t_end` be calculated instead. As a sinusoidal profile is periodic and may reach the same value multiple times, users will provide `t_end` directly and we will calculate and save the value of `X_end` in the outer constructor.

An implementation of our condition profile's struct may therefore look something like:

```@example condition_profiles
using Kinetica

mutable struct SinusoidDirectProfile{uType, tType} <: Kinetica.AbstractDirectProfile
    f::Function
    A::uType
    freq::uType
    X_start::uType
    X_end::uType
    t_end::tType
    tstops::Vector{tType}
    sol
end
nothing # hide
```

We've used a [parametric type](https://docs.julialang.org/en/v1/manual/types/#Parametric-Types) here to ensure that all values that may affect the value of our condition are using consistent numeric types.

With this definition complete, we can now define our profile's condition function and outer constructor. The condition function will take a time `t` and an instance of the `SinusoidDirectProfile`. The outer constructor will take a minimal set of parameters needed to populate and calculate the entire condition profile, and create an instance of the profile with the condition function bound to the `f` field:

```@example condition_profiles
function _f_SinusoidDirectProfile(t, profile::SinusoidDirectProfile)
    return typeof(profile.X_start)(
        ((t <= 0.0) * profile.X_start) +
        ((t > 0.0 && t <= profile.t_end) * (profile.X_start + profile.A*sin(2*pi*profile.freq*t))) + 
        ((t > profile.t_end) * profile.X_end)
    )
end

function SinusoidDirectProfile(;
    A::uType,
    freq::uType,
    X_start::uType,
    t_end::tType
) where {uType <: AbstractFloat, tType <: AbstractFloat}

    X_end = X_start + A*sin(2*pi*freq*t_end)
    tstops = [t_end]

    return SinusoidDirectProfile(_f_SinusoidDirectProfile, A, freq, X_start, X_end, t_end, tstops, nothing)
end
nothing # hide
```

How we've defined `_f_SinusoidDirectProfile(t, profile)` here may look a little strange. This is actually a way of writing if-else statements that is compatible with [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl), and is equivalent to writing
```julia
function _f_SinusoidDirectProfile(t, profile::SinusoidDirectProfile)
    if t <= 0.0
        return profile.X_start
    elseif t > 0.0 && t <= profile.t_end
        return profile.X_start + profile.A*sin(2*pi*profile.freq*t)
    else
        return profile.X_end
    end
end
```
However, Symbolics isn't compatible with this form of if-else statement so instead we use a bit of boolean algebra to embed the same meaning into our condition function. We've therefore assumed that simulation time will start at `t=0.0`, before which time the condition will be `X_start`. After `t_end`, the condition profile will stop varying and be set to the value of the sinusoid at this time. In between these two times, we take the instantaneous value of the sinusoid instead.

!!! note "Naming Condition Functions"
    Condition functions do not have to take any particular names, although Kinetica generally follows the pattern set above. Starting with an underscore indicates that the function is not intended for general use, the `f` (or `grad` in the case of gradient functions) specifies the type of variable condition profile that it attaches to, and the name of the condition profile indicates that it should only be used in the context of that particular profile.

Finally, we need to define a method of [`Kinetica.create_discrete_tstops!`](@ref) for our new profile. This function creates time stopping points for discrete rate update simulations. Some profiles contain periods of time where their condition is stationary, during which rate constant updates are unnecessary, so their method of this function can avoid creating time points during these periods. However, since the gradient of a sinusoid is only ever momentarily stationary, this optimisation is not necessary. Our implementation of this method is therefore as follows:

```@example condition_profiles
function Kinetica.create_discrete_tstops!(profile::SinusoidDirectProfile, ts_update::AbstractFloat)
    if ts_update > profile.t_end 
        throw(ArgumentError("Error defining tstops, `ts_update` is too large.")) 
    end
    profile.tstops = create_savepoints(0.0, profile.t_end, ts_update)
    return
end
nothing # hide
```

This is simply checking that the rate update timestep `ts_update` is not larger than the total duration of our condition profile, then setting the `tstops` of our profile to an evenly spaced array of time points between `t=0.0` and `t=t_end`. [`create_savepoints`](@ref) is just a Kinetica utility function for creating evenly spaced ranges of values while avoiding floating point rounding errors.

!!! warning "Adding methods to Kinetica functions"
    Note that when we extended [`Kinetica.create_discrete_tstops!`](@ref) above, we specifically prepended the function name with `Kinetica.`. The extra name of the module here is **required** as we want to add a method to Kinetica's existing function, not define our own in the `Main` module!

With this all implemented, our directly variable condition profile is now ready for use! Let's start by constructing a [`ConditionSet`](@ref) for a discrete rate update simulation with a sinusoidally-varying temperature:

```@example condition_profiles
conditions = ConditionSet(Dict(
    :T => SinusoidDirectProfile(;
        A = 75.0,
        freq = 0.25,
        X_start = 900.0,
        t_end = 20.0
    )),
    ts_update=tconvert(0.1, "ms", "s")
)
nothing # hide
```

Because we've made `SinusoidDirectProfile` a subtype of `Kinetica.AbstractDirectProfile`, Kinetica's [`ConditionSet`](@ref) has done a bunch of work in the background and set this profile up for use within Symbolics.jl by [registering](https://symbolics.juliasymbolics.org/stable/manual/functions/) the generated `SinusoidDirectProfile.f` function (which is really just a reference to `_f_SinusoidDirectProfile`) into Symbolics' computation graph. This means it's ready for use within kinetic simulations, just like that!

```@example condition_profiles
pars.tspan = (0.0, get_t_final(conditions)) # hide
# Load in existing CRN...
# pars = ODESimulationParams(...)
# calc = PrecalculatedArrheniusCalculator(...)

solvemethod = VariableODESolve(pars, conditions, calc)
res = solve_network(solvemethod, sd, rd)
nothing # hide
```

```@example condition_profiles
using Plots

p1 = plot(res)
p2 = conditionsplot(res, :T)
plot(p1, p2, layout=(2, 1))
savefig("../assets/tutorials/condition_profiles/direct.svg"); nothing # hide
```

![](../assets/tutorials/condition_profiles/direct.svg)

As you can see, the condition profile has been seamlessly integrated into Kinetica. Within [`solve_network`](@ref), a new method of [`Kinetica.solve_variable_condition!(::Kinetica.AbstractDirectProfile, ::ODESimulationParams)`](@ref) was created and run for our condition profile, which generated an interpolable `DiffEqArray` that is now available in the solved profile's `sol` field:

```@example condition_profiles
profile = get_profile(res.conditions, :T)
profile.sol
```

## Gradient-Variable Profile Implementation

Implementing a gradient-variable condition profile is a very similar procedure to the above directly variable profile - we'll still define a new struct that needs some specific fields to be compatible with the rest of Kinetica, we'll still make an outer constructor for this struct that creates a specialised function which will be automatically registered within Symbolics.jl, and we'll need a new method for [`Kinetica.create_discrete_tstops!`](@ref).

Unlike the directly variable profile, we'll need to implement a function for our sinusoid's *gradient* with respect to time. In the case of this profile, a gradient-based definition is not required since we already have an exact function for mapping our sinusiodal condition to simulation time. However, it turns out that the gradient-based representation actually helps us include an additional parameter to describe the phase of the sinusoid ``\varphi``:
```math
\frac{dX}{dt}=2\pi fA \cos \left( 2\pi f t + φ \right)
```
This would've been inconvenient to include into our directly variable form, as the user's input value of `X_start` would've needed to be changed within the outer constructor to take the starting phase into account. However, a gradient-based representation lends itself naturally to this addition.

The struct for this condition profile, which we'll call `SinusoidGradientProfile`, requires the same fields as its directly variable counterpart. The exception to this is the `f` field, which previously held the actual time-dependent condition function. In gradient-based profiles, this is replaced with a field called `grad`, which holds the gradient function. We also need to subtype `Kinetica.AbstractGradientProfile` this time, and we'll add a field for our starting phase ``\varphi``. Our new struct might therefore look something like this:

```@example condition_profiles
mutable struct SinusoidGradientProfile{uType, tType} <: Kinetica.AbstractGradientProfile
    grad::Function
    A::uType
    freq::uType
    φ::uType
    X_start::uType
    t_end::tType
    tstops::Vector{tType}
    sol
end
nothing # hide
```

The condition function and outer constructor for this profile will also be similar to their directly variable counterparts, again with an extra parameter `φ` added:

```@example condition_profiles
function _grad_SinusoidGradientProfile(t, profile::SinusoidGradientProfile)
    return typeof(profile.X_start)(
        ((t <= 0.0) * 0.0) +
        ((t > 0.0 && t <= profile.t_end) * (2*pi*profile.freq*profile.A*cos(2*pi*profile.freq*t + profile.φ))) + 
        ((t > profile.t_end) * 0.0)
    )
end

function SinusoidGradientProfile(;
    A::uType,
    freq::uType,
    φ::uType,
    X_start::uType,
    t_end::tType
) where {uType <: AbstractFloat, tType <: AbstractFloat}

    tstops = [t_end]
    return SinusoidGradientProfile(_grad_SinusoidGradientProfile, A, freq, φ, X_start, t_end, tstops, nothing)
end
nothing # hide
```

The gradient function we've defined here describes a function that is stationary (gradient of zero) before `t=0.0` and after `t=t_end`, but that varies sinusoidally in between. This sinusoidal variation now also takes into account the starting phase ``\varphi``.

The implementation of [`Kinetica.create_discrete_tstops!`](@ref) can be the same as it was for the directly variable profile, as we are describing the same underlying function:

```@example condition_profiles
function Kinetica.create_discrete_tstops!(profile::SinusoidGradientProfile, ts_update::AbstractFloat)
    if ts_update > profile.t_end 
        throw(ArgumentError("Error defining tstops, `ts_update` is too large.")) 
    end
    profile.tstops = create_savepoints(0.0, profile.t_end, ts_update)
    return
end
nothing # hide
```

And with that, we're done implementing another condition profile! Let's see it in action:

```@example condition_profiles
conditions = ConditionSet(Dict(
    :T => SinusoidGradientProfile(;
        A = 75.0,
        freq = 0.25,
        φ = pi/2,
        X_start = 975.0,
        t_end = 20.0
    )),
    ts_update=tconvert(0.1, "ms", "s")
)

solvemethod = VariableODESolve(pars, conditions, calc)
res = solve_network(solvemethod, sd, rd)
nothing # hide
```

This time, because we're using a gradient-based profile, when [`solve_network`](@ref) is called Kinetica goes away and assembles a ModelingToolkit.jl `ODESystem` and solves it with respect to time, just like it does when performing the kinetic simulation of the CRN immediately afterwards. As such, the `sol` field of our profile now contains a full `ODESolution`:

```@example condition_profiles
profile = get_profile(res.conditions, :T)
profile.sol
```

We can also look at the simulation results from using our gradient-based profile. Notice how both the concentrations and the temperature have changed since we introduced the phase parameter - now temperature starts from the top of the sinusiod rather than the middle:

```@example condition_profiles
p1 = plot(res)
p2 = conditionsplot(res, :T)
plot(p1, p2, layout=(2, 1))
savefig("../assets/tutorials/condition_profiles/gradient.svg"); nothing # hide
```

![](../assets/tutorials/condition_profiles/gradient.svg)

