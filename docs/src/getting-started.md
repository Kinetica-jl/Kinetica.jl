# Getting Started

Here we will demonstrate how to perform a simple CRN exploration using Kinetica. Simulations within Kinetica are designed to be run as a script, consisting of the following elements:

1) Parameter blocks,
2) Simulation,
3) Analysis.

We will begin by going through a basic example of each of these elements. Further details can be found in the Tutorials section of this documentation. 

If you are attempting to recreate this tutorial for yourself, it relies on a few input files that are separate to any of the main code repositories. These files can be found in the `examples` directory of this documentation's repository, which can be accessed by cloning this repository:

```bash
git clone https://github.com/Kinetica-jl/KineticaDocs.jl.git
cd KineticaDocs.jl/examples
```

To start with the tutorial, load the main Kinetica.jl package:

```@example getting_started
using Kinetica
```

We'll also set a seed for Julia's random number generation, to ensure that this tutorial is fully reproducible. This can be ignored in regular use, but is useful here if you want to compare your results to those obtained here.

```@example getting_started
using Random
Random.seed!(12345)
nothing # hide
```

## Parameter Blocks

The initial part of Kinetica simulation scripts usually consist of 3-4 blocks of parameters, depending on whether or not the simulation consists of only a kinetic calculation on a pre-existing CRN or if it also requires a CRN exploration before such a calculation can take place.

### Simulation Conditions

The first parameter block should be a [`ConditionSet`](@ref), which defines the experimental conditions a kinetic calculation should take place under. This block is usually defined first as later parameters often depend on it.

[`ConditionSet`](@ref) blocks allow for flexible definitions of any number of arbitrary conditions, each of which can be static (constant) or variable (time-dependent). For this simple CRN exploration, we will define a [`ConditionSet`](@ref) that specifies a linear temperature increase, from 300 K at time `t = 0.0` to 1000 K at a rate of 50 K/s. Kinetica comes with a library of variable condition profiles, allowing this [`ConditionSet`](@ref) to be simply defined as follows:

```@example getting_started
conditions = ConditionSet(Dict(
    :T => LinearGradientProfile(;
        rate = 50.0,
        X_start = 500.0,
        X_end = 1200.0
    )
))
```

This generates a [`LinearDirectProfile`](@ref) for the linear temperature increase we are interested in, and binds it to the symbol `:T` for use within kinetic calculators. Further information on the types of condition profiles implemented in Kinetica can be found in the tutorial on [Arbitrary Simulation Conditions](@ref).

We will be able to visualise this condition profile shortly, but we must first define a parameter block of ODE solution parameters.

### Kinetic Simulation Parameters

The [`ODESimulationParams`](@ref) block defines all of the parameters needed when converting a CRN to a system of ODEs and integrating it in time. This includes parameters such as the simulation timespan, initial concentrations of reactants and the ODE solver being used.

For the purposes of this tutorial, we will construct the following [`ODESimulationParams`](@ref) block:

```@example getting_started
using OrdinaryDiffEq
using Sundials

pars = ODESimulationParams(
    tspan = (0.0, get_t_final(conditions)),
    u0 = Dict("C" => 1.0),
    solver = CVODE_BDF(; linear_solver=:KLU)
)
```

This block details the three essential parameters for any simulation in Kinetica:

1) `tspan`: Simulation timespan (in seconds, unless otherwise specified within kinetic calculator). Must be a tuple of `(start_time, end_time)`. Here, we fetch the end time directly from the [`ConditionSet`](@ref) defined above using `get_t_final(conditions)`, which calculates the time at which all defined conditions have reached their final state.
2) `u0`: Dictionary of initial concentrations. Here we define that kinetic simulations should start with 1.0 mol dm``^{-3}`` of methane (`C` in SMILES notation), and no other reactants.
3) `solver`: ODE solution algorithm to use, from those available in [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). Here we load in the `CVODE_BDF` solver from the [Sundials.jl](https://github.com/SciML/Sundials.jl) library, as this has been the best general-purpose solver in our testing. The KLU linear solver allows this ODE solver to work with sparse arrays, which are enabled by default.

Further details of the parameters available in this block can be found in the [ODE Solution](@ref) section.

### CRN Exploration Parameters

In order to perform a CRN exploration, an exploration method must be chosen. Kinetica currently provides two: [`DirectExplore`](@ref) and [`IterativeExplore`](@ref). For now we will use the former, as it is simpler.

[`DirectExplore`](@ref) explores all chemical reactions within a given radius of the starting reactants, irrespective of whether or not they will occur within a later kinetic simulation under the selected environmental conditions. This method is best suited to small CRNs under kinetically slow conditions, where few reactions are possible and complete sampling of all available reactions is easy.

Reactions are sampled using [CDE](https://github.com/HabershonLab/cde), an external code for graph-droven sampling of reactions that is included with Kinetica.jl. CDE has its own parameters which must be set, but these are usually very similar for most CRN explorations. As such, a directory of template inputs must be provided for CDE to function. We provide such a template directory with this documentation, which can be used to run this tutorial. If you have cloned this documentation's repository, as suggested at the start of this tutorial, these files are in `KineticaDocs.jl/examples/cde_template`.

Once this is done, the exploration parameters can be set up as follows:

```@example getting_started
crn_dir = "./my_CRN"

exploremethod = DirectExplore(
    rdir_head = crn_dir,
    reac_smiles = ["C"],
    rxn_convergence_threshold = 5,
    cde = CDE(
        template_dir = "../../examples/cde_template",
        radius = 5,
        sampling_seed = 1
    )
)
```

The parameters defined in this block are as follows:

1) `rdir_head` : This is the directory in which the raw CRN exploration files will be stored. We have defined `crn_dir` above because it will be useful to reference later. 

!!! warning "I/O intensity"
    Many CDE runs can be quite I/O-heavy, so it is useful to place this directory on a fast, directly-attached drive. Running the CRN exploration out of a network-attached spinning hard drive will probably be a bad time.

2) `reac_smiles`: Array of SMILES strings representing reactants in the initial reactant system. These will be fed to CDE to generate reactions. In this case, we will generate all reactions that can occur between two methane molecules.

!!! note "The Iterative method works differently!"
    Watch out for the meaning of these parameters changing slightly between [`DirectExplore`](@ref) and [`IterativeExplore`](@ref). They can serve the same purpose, but methodologically do different things under the hood. See the tutorial on [Iterative CRN Exploration](@ref) for more information.

3) `rxn_convergence_threshold`: How many CDE iterations to perform with no new reactions discovered before considering the CRN converged.

4) `cde`: CDE parameter block using the example input templates (change this path to the directory you saved yours into). This defines a reactive radius of 5, meaning that all reactions within 5 steps of the reactant system will be generated. The `sampling_seed` parameter is usually not set, but we will set it here to make this tutorial reproducible. More information on this sub-block is also given in the [Iterative CRN Exploration](@ref) tutorial.

### Kinetic Calculator

The final parameter block to set up is the simulation's kinetic calculator. This is an object capable of calculating rate constants for every reaction in a CRN under a given set of conditions. Kinetic calculators are one of the core points of modularity in Kinetica, with extension packages like KineticaKPM.jl extending the functionality of the main code by adding in more calculators for different requirements.

The main Kinetica.jl package only includes one kinetic calculator, the [`PrecalculatedArrheniusCalculator`](@ref). This calculator calculates temperature-dependent rate constants using the Arrhenius equation:

```math
k = Ae^{-\dfrac{E_a}{RT}}
```

The [`PrecalculatedArrheniusCalculator`](@ref) requires vectors of Arrhenius prefactors `A` and activation energies `Ea` for all reactions before the code is executed, and as such is typically only used when a CRN has been generated and these values have been determined outside Kinetica. However, for the purposes of this tutorial (where the random seed for CDE's reaction generation has been set, see above), we know the reactions that are going to be generated in advance, so we have provided approximate values to input into this calculator. These values are in `KineticaDocs.jl/examples/getting_started/arrhenius_params.bson`, and can be loaded in with:

```@example getting_started
using BSON
calc_pars = BSON.load("../../examples/getting_started/arrhenius_params.bson")
```

again, replacing this path with the equivalent path on your computer. The calculator for this CRN can now be constructed:

```@example getting_started
calc = PrecalculatedArrheniusCalculator(calc_pars[:Ea], calc_pars[:A]; k_max=1e12)
```

For more advanced calculators that allow on-the-fly calculation of rate constants as reactions are generated, see the tutorial on [Kinetic Calculators](@ref).

## Simulation

Now that all of the parameter blocks have been constructed, we can start generating and simulating CRNs!

### Visualising Condition Profiles

First, now that we have both a [`ConditionSet`](@ref) and a [`ODESimulationParams`](@ref) block, we can have a look at the variable temperature profile that we created. To begin, we need to construct a DifferentialEquations `ODESolution` for this profile. Kinetica provides a shortcut for this, with a one-liner that constructs `ODEProblem`s and solves them for every variable condition in a [`ConditionSet`](@ref):

```@example getting_started
solve_variable_conditions!(conditions, pars)
```

!!! note "Not neccessary in regular use!"
    We are manually triggering the solution of the [`ConditionSet`](@ref) here to visualise how the constructed temperature profile will look in the final simulation. However, calling [`solve_variable_conditions!`](@ref) as above is not neccessary in all Kinetica scripts. This function is normally run automatically before integrating the CRN ODEs, and does not typically need to be called like this.

With this in place, we can now plot our temperature profile's `ODESolution` to see how the simulation temperature is going to vary with time:

```@example getting_started
mkpath("assets/getting-started") # hide
using Plots
plot(get_profile(conditions, :T).sol)
savefig("assets/getting-started/Tprofile.svg"); nothing # hide
```

![](assets/getting-started/Tprofile.svg)

As expected, the temperature profile we constructed goes from 300 K to 1000 K at a rate of 50 K/s. Provided we are happy with this, we can continue to perform a CRN exploration and a kinetic simulation.

### Running the Simulation

CRN exploration and kinetic simulation are wrapped under a single function call: [`explore_network`](@ref). This takes the exploration parameters, the simulation parameters and calculator (which we will wrap under a single variable), and optionally a directory to save the finished results to. These results will additionally be returned by the function. To run the entire simulation, we simply call the following:

```@example getting_started
solvemethod = VariableODESolve(pars, conditions, calc)
mkdir("./my_CRN_out") # hide
res = explore_network(exploremethod, solvemethod; savedir="./my_CRN_out")
nothing # hide
```

!!! note "Checkpoints and Restarts"
    In the event that a Kinetica CRN exploration fails before completion, all is not lost! As long as you still have the original head directory that the CRN was being explored in (`crn_dir` in this tutorial), both of Kinetica's exploration algorithms will detect where you left off and restart from there seamlessly. In the iterative exploration algorithm, Kinetica also creates checkpoint files with current CRN state and kinetic simulation results, which are saved under the `savedir` keyword argument specified above along with the final simulation output.

The resulting [`ODESolveOutput`](@ref) object contains the explored CRN, the kinetic simulation results, and all of the parameters and conditions that went into it. These can be easily and efficiently saved and loaded as needed (see [Saving & Loading](@ref)).

## Analysis

The [`ODESolveOutput`](@ref) object has plot recipes defined for easy plotting using [Plots.jl](https://github.com/JuliaPlots/Plots.jl), allowing complex figures and statistics to be shown with a single line of code. To begin, we can look at how the concentrations of all of the species found during CRN generation vary over time under the variable temperature profile we specified. This is achieved by simply running

```@example getting_started
plot(res)
savefig("assets/getting-started/kinetics_plot.svg"); nothing # hide
```

![](assets/getting-started/kinetics_plot.svg)

This CRN presents some interesting results! Accoding to the kinetics enforced by the calculator we have used, methane will not start breaking down into any of the free radical species discovered within this CRN until ``t \approx 5 \text{ s}``. We could check the temperature that this occurs at by referencing the temperature profile we plotted above but this is also accessible from `res` by running

```@example getting_started
conditionsplot(res, :T)
savefig("assets/getting-started/Tprofile_2.svg"); nothing # hide
```

![](assets/getting-started/Tprofile_2.svg)

Notice how, because we used the `conditionsplot()` function, the correct symbolically-indexed condition profile was obtained and axis labels were set accordingly.

If we want to obtain the numeric value of the temperature at ``t \approx 5 \text{ s}``, we can interpolate this directly:

```@example getting_started
Tprofile = get_profile(res.conditions, :T)
Tprofile.sol(5.0)
```

If we had multiple variable condition profiles, the above could also be done for them by simply passing in their bound symbols.

Similarly, we can also obtain the concentrations of all of the species in the reaction mixture at any time in the simulation through interpolation:

```@example getting_started
res.sol(5.0)
```

It is also very useful to analyse the final concentrations of species at the end of a kinetic simulation. Kinetica defines another plot recipe for this:

```@example getting_started
finalconcplot(res)
savefig("assets/getting-started/concs_plot.svg"); nothing # hide
```

![](assets/getting-started/concs_plot.svg)

## Next Steps

Now you're familiar with the basics of network exploration and kinetic simulation within Kinetica, you can learn more by looking through the other tutorials! 
