# Calculator Interface

For an introduction to Kinetica's kinetic calculators, which allow rate constant calculation under arbitrary experimental conditions of interest, see the tutorial on [Kinetic Calculators](@ref). This page details the underlying implementation of these calculators in order to facilitate extension and development of new calculators.

All kinetic calculators implemented in Kinetica and its modules follow the same basic interface. They should consist of:

* A (mutable) `struct` which holds all of the data that the calculator needs to calculate rate constants. This must be a subtype of `Kinetica.AbstractKineticCalculator`. Mutability is not a prerequisite, but is suggested as many calculators are designed to be updated further down the calculation chain.
* An outer constructor which simplifies construction of the main struct, potentially offering optional arguments such as the ability to change the time unit the calculator works in. While not strictly necessary, this is recommended.
* A method of the [`setup_network!`](@ref) function which can be called with the calculator being implemented. This handled compatibility checking with the current CRN, as well as any additional calculation that needs to take place on a per-reaction basis (and therefore cannot be done during calculator construction).
* A method of [`splice!(calc::cType, rids::Vector{Int}) where cType<:Kinetica.AbstractKineticCalculator`](@ref) which can be called with the calculator being implemented. This enables data about unnecessary reactions to be deleted, mirroring the funcitonality of [`splice!(::RxData, ::Vector{Int})`](@ref) for deletion of reactions from the current CRN.
* A method of [`allows_continuous`](@ref) for the calculator being implemented, which specifies if the calculator is allowed to be used in continuous rate update simulations (see [ODE Solution](@ref)).
* A method of [`has_conditions`](@ref) for the calculator being implemented, which checks if the experimental conditions in the current kinetic simulation are compatible with the calculator.
* Functors of the calculator being implemented which calculate the rate constant of every reaction in the current CRN at a given set of experimental conditions. These must accept the experimental conditions of interest as keyword arguments.

As long as these are implemented, a new calculator will slot in seamlessly to Kinetica's arbitrary variable condition and CRN simulation frameworks.

## Implementation Example

We will demonstrate the implementation of a new kinetic calculator with a somewhat contrived example that should nevertheless demonstrate all of the above points. This calculator will return rate constants that:

* Are dependent on temperature `:T` and pressure `:P` as experimental conditions,
* Can be calculated with or without a maximum rate constant `k_max`,
* Take a base rate constant `k_base` and multiply it by the ratio of products to reactants in a given reaction, under the formula
```math
k_r = \frac{n_p}{n_r}k_{\text{base}} \exp\left[-\frac{P}{10T}\right]
```
where ``n_p`` is the number of product molecules in a given reaction, and ``n_r`` is the number of reactant molecules.

* Assuming the above equation yields rate constants in time units of ``s^{-1}``, can be modified to work under other time units.

This rate constant contains no physical meaning, but will illustrate all of the requirements of a kinetic calculator very well.

### Calculator Struct

First, we must create the main calculator struct, which will hold user-facing parameters such as `k_base` and `k_max` that are set at calculator construction and are not intended to be modified much down the line, as well as internal parameters such as `n_reacs` and `n_prods` which are reaction-specific and will be calculated later. Assuming we're naming this calculator `MyNewCalculator`, this can be done as follows:

```@example calculator_interface
using Kinetica

mutable struct MyNewCalculator{kmType, uType, tType} <: Kinetica.AbstractKineticCalculator
    k_base::uType
    n_reacs::Vector{Int}
    n_prods::Vector{Int}
    k_max::kmType
    t_unit::String
    t_mult::tType
end
```

Here we've defined a [Parametric Type](https://docs.julialang.org/en/v1/manual/types/#Parametric-Types), which allows us to both specialise to the numeric types of species concentrations `uType` (assuming we represent `k_base` with the same numeric type as our species concentrations) and of simulation time `tType`, as well as allowing us to use multiple dispatch based on whether or not we have a `k_max` using `kmType` down the line. These parameters are not necessary to create a functioning calculator, but we include them here for consistency with how most of Kinetica's calculators are curently implemented.

### Outer Constructor

While it's nice to keep the information for both what our time unit is `t_unit` and what it does to our rates `t_mult` for reference, one quantity can be derived from the other and therefore both do not need to be provided to construct our calculator. We may also want to establish some sensible defaults for both `k_max` and `t_unit` to avoid having to input them every time we need to construct a calculator. Finally, we also don't need to provide values to CRN-specific fields such as `n_reacs` and `n_prods` at construction time. These features can be achieved by implementing an outer constructor:

```@example calculator_interface
function MyNewCalculator(k_base::uType; k_max::Union{uType, Nothing}=nothing, 
                         t_unit="s") where {uType <: AbstractFloat}
    
    t_mult = tconvert(t_unit, "s")
    return MyNewCalculator(k_base, Int[], Int[], k_max, t_unit, t_mult)
end
nothing # hide
```

With this, we've now created a default state so that when a user calls `calc = MyNewCalculator(k_base)`, the calculator creates empty arrays for `n_reacs` and `n_prods`, and automatically assumes there is no upper rate constant limit and we're working in rate units where time is measured in seconds. This instance of the calculator will be of type `MyNewCalculator{Nothing, uType, tType}`, allowing us to dispatch on a rate constant equation which excludes `k_max`. We've also enforced that the type of `k_max` (when provided) must be the same as that of `k_base` by making this constructor a [Parametric Method](https://docs.julialang.org/en/v1/manual/methods/#Parametric-Methods), and that this type must be a subtype of `AbstractFloat`.

### Simulation Setup

So far we've taken care of the parameters that are specific to a given instance of our calculator, but what about the parameters that are specific to a CRN being simulated? In the above outer constructor, we've left `n_reacs` and `n_prods` (which represent the number of reactant and product molecules in each reaction) as empty integer arrays. We need these values to calculate our rate constants, so we'll need to populate these arrays within [`setup_network!`](@ref).

[`setup_network!`](@ref) is a function that is always called early in a kinetic simulation under [`solve_network`](@ref). It takes the current CRN as a set of [`SpeciesData`](@ref) and [`RxData`](@ref) objects, along with the given kinetic calculator, and modifies the internal state of the calculator to work with the CRN - in our example here, it will calculate and update `n_reacs` and `n_prods`:

```@example calculator_interface
function setup_network!(sd::SpeciesData, rd::RxData, calc::MyNewCalculator)
    n_reacs = zeros(Int, rd.nr)
    n_prods = zeros(Int, rd.nr)
    # Number of reactant/product molecules is the sum of the 
    # reactant/product stoichiometries for each reaction.
    for i in 1:rd.nr
        n_reacs[i] = sum(rd.stoic_reacs[i])
        n_prods[i] = sum(rd.stoic_prods[i])
    end

    calc.n_reacs = n_reacs
    calc.n_prods = n_prods
    return
end
nothing # hide
```

!!! warning "Fixed Argument Structure"
    As [`setup_network!`](@ref) is intended to be called automatically as part of a chain of functions, it must always use the same call signature of `(sd, rd, calc)`, even when one or more of these arguments is unused, as `sd` is here. Similarly, even if no calculator setup is required, an empty method **must** be implemented here so that normal CRN simulation can proceed.

### Calculator Splicing

Sometimes, reactions are removed from CRNs as part of the kinetic simulation workflow (e.g. as part of a low-rate reaction cleanup, see `low_k_cutoff` in [ODE Solution](@ref)). Since this occurs **after** the calculator setup done in [`setup_network!`](@ref), any reactions removed from a CRN's [`RxData`](@ref) through [`splice!(::RxData, ::Vector{Int})`](@ref) must similarly be removed from any reaction-specific fields in the calculator (in our case, `n_reacs` and `n_prods`). This is done through a custom method of [`splice!(calc::cType, rids::Vector{Int}) where cType<:Kinetica.AbstractKineticCalculator`](@ref):

```@example calculator_interface
function Base.splice!(calc::MyNewCalculator, rids::Vector{Int})
    # Splice the underlying reaction-specific arrays.
    splice!(calc.n_reacs, rids)
    splice!(calc.n_prods, rids)
    return
end
nothing # hide
```

!!! warning "Required Implementation"
    As with [`setup_network!`](@ref), a method of [`splice!(calc::cType, rids::Vector{Int}) where cType<:Kinetica.AbstractKineticCalculator`](@ref) **must** be implemented for every kinetic calculator, as it can be called as part of [`solve_network`](@ref). 

### Continuous Rate Tagging

Some calculators may have rate constant calculation functions that contain calls to external programs or are otherwise incompatible with representation as a [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)-based expression through [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl). These calculators are not compatible with the continuous rate update formalism and must be run through discrete rate update simulations instead (see [ODE Solution](@ref) for details).

Kinetica prevents these calculators from being passed to an incompatible solution algorithm by checking the calculator's [`allows_continuous`](@ref) method. This simply returns a `Bool` stating whether the calculator should be allowed within continuous formalism simulaitons. In our example here, the rate constant function that we are going to implement can be represented in symbolic algebra by Symbolics, so we simply add:

```@example calculator_interface
allows_continuous(::MyNewCalculator) = true
nothing # hide
```

!!! note "Which rate expressions are allowed?"
    If you're new to Symbolics or ModelingToolkit, it may be unclear which rate expressions are or are not compatible with the continuous formalism. In some cases this may require some trial and error, but in most, if you can write out the algebraic form of the rate expression by hand, it can be made to be compatible. If you need to use loops, this won't work (vectorised operations are okay though). Don't worry if your rate expression isn't compatible - the discrete rate update formalism is usually much quicker to compile, and sometimes it's quicker to solve too!

### Argument Pre-Checking

It's often useful to check far in advance if a kientic calculator is actually capable of working with the experimental conditions you're throwing at it. This is performed automatically by a method of [`has_conditions`](@ref) when the calculator and the [`ConditionSet`](@ref) first meet, in the inner constructor of [`StaticODESolve`](@ref)/[`VariableODESolve`](@ref).

We must therefore implement such a method for our calculator, which lets Kinetica know which symbolic conditions are acceptable:

```@example calculator_interface
function has_conditions(::MyNewCalculator, symbols::Vector{Symbol})
    return all([sym in [:T, :P] for sym in symbols])
end
nothing # hide
```

This way, if we were to provide a `:V` condition in our [`ConditionSet`](@ref) which our calculator wouldn't know what to do with, this would be caught early instead of causing problems deep into a difficult to diagnose [`solve_network`](@ref) call.

### Rate Constant Calculation

Finally, with the rest of the calculator out of the way, we can write our actual rate expressions. These are written as [functors](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects) of our main calculator struct. Here, we will define two functors and allow Julia to dispatch to the correct one based on the type of calculator it receives. If the calculator's `kmType` is the same as its `uType`, we know that a value was provided for `k_max` and we need to dispatch on a modified version of the rate expression. Otherwise, if `kmType` is `Nothing`, we know there's no `k_max` and we can dispatch on the normal rate expression:

```@example calculator_interface
# Dispatching without k_max awareness.
function (calc::MyNewCalculator{Nothing, uType, tType})(; T::Number, P::Number) where {uType, tType}
    k_r = (calc.n_prods ./ calc.n_reacs) .* calc.k_base * exp(-P/(10*T)) * calc.t_mult
    return k_r
end

# Dispatching with k_max awareness.
function (calc::MyNewCalculator{uType, uType, tType})(; T::Number, P::Number) where {uType, tType}
    k_r = (calc.n_prods ./ calc.n_reacs) .* calc.k_base * exp(-P/(10*T)) * calc.t_mult
    return 1.0 ./ ((1.0 / calc.k_max) .+ (1.0 ./ k_r))
end
nothing # hide
```

Vitally important here are:
* Experimental conditions must be passed in as **keyword arguments**, as Kinetica splats an array of `Pair`s of condition `Symbol`s and their values (static parameters or variable expressions) into the calculator internally to enable arbirtary conditions to be passed to arbitrary calculators. This is accomplished by indicating that there are no positional arguments with a semicolon before listing the conditions.
* Calculator functors return rate constants for **all** reactions in a CRN simultaneously, so rate constant expressions typically make use of vectorised operations to make neat one-liners like those above.

### Checking Our Work

That's it! Now we have everything implemented for a working (albeit non-physical and mostly useless) kinetic calculator. We can make sure it works by running all the above functions on a small CRN. We'll use the CRN that was created for [Getting Started](@ref):

```@example calculator_interface
res = load_output("../my_CRN_out/direct_network_final.bson")
sd, rd = res.sd, res.rd

# Test added methods.
calc = MyNewCalculator(2e5)
println("Calculator allows continuous formalism: $(allows_continuous(calc))")
println("Calculator accepts conditions Z, P, V: $(has_conditions(calc, [:Z, :P, :V]))")
println("Calculator accepts conditions T, P: $(has_conditions(calc, [:T, :P]))")
setup_network!(sd, rd, calc)
@assert length(calc.n_reacs) == length(calc.n_prods) == rd.nr
println("Pre-splice calculator contains info on $(rd.nr) reactions")
splice!(rd, calc, collect(1:5)) # Remove first 5 reactions from rd and calc
@assert length(calc.n_reacs) == length(calc.n_prods) == rd.nr
println("Post-splice calculator contains info on $(rd.nr) reactions")

# Default calculator implementation.
calc = MyNewCalculator(2e5)
setup_network!(sd, rd, calc)
k_default = calc(; T=1000.0, P=1e5)

# Calculator with rate constants capped to 10 s-1
calc = MyNewCalculator(2e5; k_max=10.0)
setup_network!(sd, rd, calc)
k_capped = calc(; T=1000.0, P=1e5)

# Compare default and capped rate constants.
using Printf
println("\nRID  k_default  k_capped")
for i in 1:rd.nr
    println("$(rpad(string(i), 3, ' '))  $(@sprintf("%.9f", k_default[i]))  $(@sprintf("%.9f", k_capped[i]))")
end
nothing # hide
```