abstract type AbstractSimulationParams end

"""
Container for ODE-driven simulation parameters.

Contains fields for:
* Simulation timespan (must match time unit used in an attached calculator) (`tspan`)
* Initial concentrations of species, either as a Dict of certain species or a Vector of all species (`u0`)
* DifferentialEquations ODE solver (`solver`)
* Whether to use ModelingToolkit to formulate an analytical Jacobian (`do_jac=true`)
* Whether to use ModelingToolkit to formulate a sparse problem (`do_sparse=true`)
* Absolute tolerance of ODE solver (`abstol=1e-10`)
* Relative tolerance of ODE solver (`reltol=1e-8`)
* Whether to use adaptive solver tolerance (`adaptive_tols=true`)
* Whether to update solver tolerances after successful solve with adaptive tolerance (`update_tols=false`)
* Whether to break solution into chunks of size `solve_chunkstep` to avoid floating point underflow (`solve_chunks=true`)
* Global timestep at which solution should be reinitialised when `solve_chunks=true` (`solve_chunkstep=1e-3`)
* Maximum number of ODE solver iterations (`maxiters=1e5`)
* Whether to explicitly disallow negative values in the solver (`ban_negatives=false`)
* Whether to display progress bars - requires `TerminalLogger` initialisation (`progress=false`)
* Time interval to interpolate solution data on (`save_interval=nothing`)
* Cutoff below which reactions with low rate constants are removed from the network (`low_k_cutoff=:auto`)
* Whether to allow a vector `u0` to be shorter than the number of species in the network (`allow_short_u0=false`)
"""
@kwdef mutable struct ODESimulationParams{tType, uType} <: AbstractSimulationParams
    # Main parameters that dictate the entire simulation
    tspan::Tuple{tType, tType}
    u0::Union{Dict{String, uType}, Vector{uType}}

    # Solver parameters
    solver
    jac::Bool=true
    sparse::Bool=true
    abstol::uType = 1.0e-10
    reltol::uType = 1.0e-8
    adaptive_tols::Bool = true
    update_tols::Bool = false
    solve_chunks::Bool = true
    solve_chunkstep::tType = 1e-3
    maxiters::Integer = 100000
    ban_negatives::Bool = false
    progress::Bool = false

    # Optional network parameters
    save_interval::Union{tType, Nothing}=nothing
    low_k_cutoff::Union{uType, Symbol}=:auto
    allow_short_u0::Bool=false
end