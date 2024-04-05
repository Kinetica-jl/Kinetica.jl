# Kinetic Calculators

Kinetica allows for any and all rate constant expressions to be used within kinetic simulations through its modular kinetic calculator interface. This allows users and developers alike to quickly define a set of functions which are dependent on a CRN (through a [`SpeciesData`](@ref) and an [`RxData`](@ref)) and some arbitrary experimental conditions which the simulation takes place under, such as temperature and pressure, and to calculate rate constants. For more information on how kinetic calculators are implemented, see [Calculator Interface](@ref).

Kinetica features a single base calculator, the [`PrecalculatedArrheniusCalculator`](@ref). This calculator relies on having Arrhenius prefactors and activation energies precalculated for every reaction in a CRN, and serves mostly as a test calculator and an implementation example. 

Other more useful calculators are provided as modular addon packages that extend Kinetica.jl (e.g. KineticaKPM.jl). Some calculators have extensive dependencies, so this modularisation allows for picking and choosing only the calculators required for a specific project rather than having every available calculator in one bloated Julia/Python environment.

## Manually Calling Calculators

The examples below all show how calculators can be manually called to calculate rate constants. This is useful when these values are required, e.g. for analysis, but is not required during kinetic simulations. Passing a calculator into [`explore_network`](@ref) or [`solve_network`](@ref) via a [`VariableODESolve`](@ref) struct (see the ['Running the Simulation' section of Getting Started](@ref "Running the Simulation")) is enough to set up the calculator and evaluate rate constants as many times as is required by the simulation.

In the event that a calculator is required to be called manually, the process is very simple. Once the calculator has been instantiated, it should be passed to its `setup_network!` method. All kinetic calculators implement such a function, which checks the values provided to the calculator are compatible with a given [`SpeciesData`](@ref) and [`RxData`](@ref) and does any necessary pre-calculation such that the calculator can perform the minimum computation when rate constants are requested. Once setup, the calculator object can be called as a functor with its implemented experimental conditions as keyword arguments to evaluate the rate constants of all reactions in the given CRN at the provided values of the conditions. Examples of this process can be found below.

## Calculator Showcase

In all of the examples below, `sd` and `rd` refer to instances of [`SpeciesData`](@ref) and [`RxData`](@ref) respectively. These are the internal representations of species and reactions within Kinetica, see the page on [CRN Representation](@ref) for further information.

!!! note "This section is growing!"
    Kinetica currently has only a handful of kinetic calculators available. We are adding more as we need them, but if you require a specific implementation then please let us know on our [Issues page](https://github.com/Kinetica-jl/Kinetica.jl/issues). Alternatively, calculators aren't too hard to implement yourself, and custom calculators can be dropped into kinetic simulations just like the ones presented here. See [Calculator Interface](@ref) for details.

### Kinetica.jl

#### [`PrecalculatedArrheniusCalculator`](@ref)

This calculator is dependent on temperature as an experimental condition, and takes a vector of Arrhenius prefactors and a vector of activation energies, each with an entry for each reaction in a given CRN, and calculates rate constants with the Arrhenius equation. It also accepts an optional maximum rate constant `k_max` which takes over through partial diffusion control. The calculated rate constant for reaction ``i`` is therefore: 

```math
k_i = \dfrac{1}{\dfrac{1}{k_{\text{max}}} + \dfrac{1}{A_ie^{-\dfrac{E_i}{RT}}}}
```

where ``A_i`` is the ``i``th Arrhenius prefactor, ``E_i`` is the ``i``th activation energy, ``R`` is the ideal gas constant and ``T`` is the temperature, which must be passed in as a parameter.

##### Example:

```julia
# Get Arrhenius prefactors A and activation energies Ea from elsewhere...
# Length of A and Ea should match rd.nr.
calc = PrecalculatedArrheniusCalculator(Ea, A; k_max = 1e12)
setup_network!(sd, rd, calc)
k = calc(; T = 300.0)
```

### KineticaKPM.jl

The calculators in KineticaKPM all use [KineticPredictorModel](https://github.com/joegilkes/KineticPredictorModel) (KPM), a Python code for predicting activation energies using a simple neural network, as a driver for rate constant calculation using the Arrhenius equation. Arrhenius prefactors are estimated through a variety of methods depending on the calculator used.

The calculators in this package require an instance of [`KPMRun`](@ref), which acts as an interface to the KPM package and handles conversion of the current CRN into reactions which it can predict activation energies for. This is simply constructed by calling the following:

```julia
kpm = KPMRun("/path/to/kpm_model.npz")
```

where the KPM model `.npz` file should be obtained separately through training a model with the main Python package.

All calculators in this package accept an optional maximum rate constant `k_max` which takes over through partial diffusion control. They are also all capable of returning rate constants with [Measurements.jl](https://juliaphysics.github.io/Measurements.jl/stable/) uncertainties, derived from the standard deviation between activation energy predictions within an ensemble of neural networks. While these can be used manually, they are not currently supported within any of the ODE solution methods in Kinetica.

#### [`KPMBasicCalculator`](@ref)

This calculator is dependent on temperature as an experimental condition. It estimates Arrhenius prefactors as

```math
A = \dfrac{RT}{h}
```

for all reactions, where ``R`` is the ideal gas constant, ``T`` is the temperature, which must be passed as a parameter, and ``h`` is the Planck constant. The resulting rate constant for reaction ``i`` is therefore:

```math
k_i = \dfrac{1}{\dfrac{1}{k_{\text{max}}} + \dfrac{1}{\dfrac{RT}{h}e^{-\dfrac{E_i}{RT}}}}.
```

##### Example:

```julia
# Set up KPMRun object before this...
calc = KPMBasicCalculator(kpm; uncertainty = false, k_max = 1e12)
setup_network!(sd, rd, calc)
k = calc(; T = 300.0)
```

#### [`KPMCollisionCalculator`](@ref)

This calculator is dependent on temperature as an experimental condition. It estimates Arrhenius prefactors using collision theory, a hard-sphere approximation of collision frequency. All unimolecular reactions therefore require a collision partner for reactions to occur, which can be passed through the `inert_species` argument. If this is given, `setup_network!` will modify all unimolecular reactions to become bimolecular with the provided collision partners. Otherwise, an average collision partner will be calculated from the species in the CRN and rates will be calculated assuming a concentration of 1 mol dm``^{-3}`` of this 'species'.

The calculator computes two properties for each reaction: the reduced mass ``\mu`` and the collision cross section ``\sigma``. For collision partners ``A`` and ``B``, these are defined as

```math
\mu = \dfrac{m_A m_B}{m_A + m_B} \\
\sigma = \pi \left( r_A + r_B \right)^2
```

where ``m_A`` and ``m_B`` are the masses of ``A`` and ``B`` respectively and ``r_A`` and ``r_B`` are the hard sphere radii of ``A`` and ``B`` respectively.

The calculator optionally allows for specification of a 'steric factor' ``\rho``. Collision theory is known to overestimate rate constants, but there is no perfect mathematical relationship that can calculate how much this overestimate is by. The steric factors implemented here attempt to establish an empirical correction to the collision theory rate constant based on a variety of information about each species. See [`KineticaKPM.calc_steric_factors`](@ref) for further information.

The resulting rate constant for reaction ``i`` is therefore:

```math
k_i = \dfrac{1}{\dfrac{1}{k_{\text{max}}} + \dfrac{1}{\sigma_i \rho_i N_A \sqrt{\dfrac{8 k_b T}{\pi \mu_i}}e^{-\dfrac{E_i}{RT}}}}
```

where ``N_A`` is Avogadro's number, ``k_b`` is the Boltzmann constant, ``T`` is the temperature, which is passed as a parameter, ``E_i`` is the activation energy and ``R`` is the ideal gas constant.

##### Example:

```julia
# If using inert species, these need an initial concentration.
pars = ODESimulationParams(
    tspan = (0.0, 10.0)
    u0 = Dict(
        "CC" => 1.0,
        "N#N" => 1.0
    ),
    solver = ...
)
# Set up KPMRun object before this...
calc = KPMCollisionCalculator(
    kpm,
    inert_species = ["N#N"],
    steric_factor = :basic,
    uncertainty = false,
    k_max = 1e12
)
setup_network!(sd, rd, calc)
k = calc(; T = 300.0)
```