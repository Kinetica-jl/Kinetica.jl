# Implementation Details

The purpose of this page is to expand on the definition and implementation of some of the features in Kinetica that appear in other pages within this documentation. Sections here are not critical to the understanding of how to use parts of Kinetica, but may be necessary for development purpooses.

## [Chunkwise Time](@id implementation_chunkwise_time)

Chunkwise simulation time is mentioned primarily in the tutorial on [ODE Solution](@ref). It is a workaround to a floating point underflow issue that is caused by accumulating extremely small timesteps onto a relatively large global simulation time (GST), ``t_{\text{global}}``. 

To be more precise, assuming ``t_{\text{global}}`` is represented by an IEEE 64-bit floating point number (i.e. a `Float64`) where the smallest precision that is computationally representable ``\epsilon = 2^{-53} \simeq 10^{-16}``. The smallest value that can be added at any given time is therefore ``t_{\text{global}} \epsilon``. When very small timesteps are being taken and the GST is large, underflow occurs.

In practice, the impact of this problem can be seemingly minimal, with most simulations running without issue. However, closer inspection of the timesteps being taken reveals that floating point underflow can still occur in otherwise convergent simulations. We theorise that this is actually due to x86-64 processors allowing for a limited set of computations to take place on [extended precision](https://en.wikipedia.org/wiki/Extended_precision#x86_extended_precision_format) 80-bit floating point numbers where $\epsilon=2^{63}\approx10^{-18}$ when necessary, which only appear to underflow when being converted back to `Float64` in memory. However, this is slower than performing the respective 64-bit operations and is entirely implementation and platform-dependent, should not be relied upon to always work.

Enabling chunkwise time (`ODESimulationParams.solve_chunks=true`) therefore splits the overall simulation timespan (`ODESimulationParams.tspan`) into 'chunks' of local time, each of length `ODESimulationParams.solve_chunkstep`, or ``\tau_c``.

!!! warning "Uneven Chunk Splits"
    Currently, Kinetica only supports splitting simulations into chunks with the same ``\tau_c``. This means that `ODESimulationParams.tspan[2] - ODESimulationParams.tspan[1]` must be evenly divisible by `ODESimulationParams.solve_chunkstep`, such that an integer number of chunks are created.

The first simulation chunk is initialised with species concentrations ``\mathbf{c}=\mathbf{c_0}``, ``t_{\text{global}}=0.0``, local simulation time (LST) ``t_{\text{local}}=0.0`` and number of chunks ``n_c = 0``. Integration of species concentrations then proceeds with respect to ``t_{\text{local}}`` until ``t_{\text{local}}=\tau_c``, at which point the concentrations and their respective LSTs are saved to a global-time array. In this array, ``t_{\text{global}}`` is calculated as
```math
t_{\text{global}} = t_{\text{local}} + \tau_c n_c
```
and species concentrations and times are interpolated onto a grid of save times
```math
\mathbf{t_{\text{saves}}}=\left( \tau_{\text{save}}n \right)_{n=0}^{\tau_{c}/\tau_{\text{save}}}
```
where ``\tau_{\text{save}}`` is simply `ODESimulationParams.save_interval`.

The solver is then reinitialised at ``t_{\text{local}}=0.0`` with ``\mathbf{c}=\mathbf{c_{\text{final}}}``, where ``\mathbf{c_{\text{final}}}`` is the species concentrations at the end of the previous chunk, and ``n_c`` is incremented. This continues until ``t_{\text{global}}`` reaches `ODESimulationParams.tspan[2]`. By replacing ``t_{\text{global}}`` with a ``t_{\text{local}}`` that can never become very large during timestep accumulation, timesteps can be safely taken down to values of ``\tau_c \epsilon``. This allows long-timescale, high-rate simulations to proceed unhindered by floating point underflow.

## [Adaptive Solver Tolerance](@id implementation_adaptive_tolerance)

When solving `ODESystem`s, adaptive timestepping ODE solvers requires two tolerances (`abstol` and `reltol`, given within Kinetica through `ODESimulationParams.abstol` and `ODESimulationParams.reltol` respectively) which control the numeric error allowed within the solver. A balance must be struck with these tolerances, as setting very low values increases solution accuracy, at the cost of taking many more small timesteps in order to achieve this.

By default, these tolerances are set at the start of a simulation and remain constant throughout. This can be a problem if the solver runs into an area of the simulation where an accurate solution cannot be obtained at the current tolerances, as the solver will throw an error and crash. To avoid this, Kinetica implements an adaptive-tolerance solution method based around DifferentialEquations.jl's [integrator interface](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/). This method is enabled with the `ODESimulationParams.adaptive_tols` parameter.

This method takes a pre-assembled integrator, and attempts to solve the underlying `ODEProblem` at the initially requested tolerances. If this is not possible, `abstol` and `reltol` are decreased by a factor of 10, reducing numerical error at the cost of a more expensive simulation. The tolerances are decreased until one or both fall below the precision of the numeric type used within the simulation, at which point an error will be returned as the simulation cannot be solved. An error is also returned if more than five tolerance changes have been attempted, to avoid needlessly running expensive, impossible simulations.

If the solver tolerances need to be modified during adaptive tolerance solution and `ODESimulationParams.update_tols=true`, the changes to the tolerances will be echoed back to the [`ODESimulationParams`](@ref) object returned within the simulation's [`ODESolveOutput'](@ref), which can be useful if similar simulations need to be rerun with the same parameters, as the solver tolerances will stay at workable values.

!!! note "Adaptive tolerance in chunkwise simulations"
    The behaviour enabled with `ODESimulationParams.update_tols=true` is usually not desired when running simulations with chunkwise time. Decreased solver tolerances are often only required in localised areas of time within kinetic simulations where reactions are proceeding exceptionally quickly and there are very large changes in species concentrations. As such, decreased tolerances may only be necessary within a few simulation chunks, while the rest of the simulation can get by with a higher tolerance (and therefore a faster solve time within these chunks). By not updating the tolerances within [`ODESimulationParams`](@ref), each simulation chunk starts at the default tolerances and only reduces them if necessary.

## [Removing Low-Rate Reactions](@id implementation_low_rate)