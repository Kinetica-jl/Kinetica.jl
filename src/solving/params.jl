abstract type AbstractSimulationParams end

mutable struct ODESimulationParams{tType, uType} <: AbstractSimulationParams
    # Main parameters that dictate the entire simulation
    tspan::Tuple{tType, tType}
    u0::Union{Dict{String, uType}, Vector{uType}}

    # Solver parameters
    solver
    jac::Bool
    sparse::Bool
    abstol::uType
    reltol::uType
    adaptive_tols::Bool
    update_tols::Bool
    solve_chunks::Bool
    solve_chunkstep::tType
    maxiters::Integer
    ban_negatives::Bool
    progress::Bool

    # Optional network parameters
    save_interval::Union{tType, Nothing}
    low_k_cutoff::Union{uType, Symbol}
    low_k_maxconc::uType
    allow_short_u0::Bool
end

"""
Keyword-defined container for ODE-driven simulation parameters.

Catches common errors in the simulation parameters early, before
any expensive calculatons are run.

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
* Maximum species concentrationto multiply maximum reaction rates by un low rate cutoff (`low_k_maxconc=2.0`)
* Whether to allow a vector `u0` to be shorter than the number of species in the network (`allow_short_u0=false`)
"""
function ODESimulationParams(;
        tspan::Tuple{tType, tType},
        u0::Union{Dict{String, uType}, Vector{uType}},
        solver,
        jac::Bool=true,
        sparse::Bool=true,
        abstol::uType=1.0e-10,
        reltol::uType=1.0e-8,
        adaptive_tols::Bool=true,
        update_tols::Bool=false,
        solve_chunks::Bool=true,
        solve_chunkstep::tType=1e-3,
        maxiters::Integer=100000,
        ban_negatives::Bool=false,
        progress::Bool=false,
        save_interval::Union{tType, Nothing}=nothing,
        low_k_cutoff::Union{uType, Symbol}=:auto,
        low_k_maxconc::uType=2.0,
        allow_short_u0::Bool=false
    ) where {tType, uType}

    # Check that the time span is valid
    if tspan[1] >= tspan[2]
        throw(ArgumentError("Invalid time span: Start = $(tspan[1]), End = $(tspan[2])"))
    end

    # Check low_k_cutoff is a valid value
    if typeof(low_k_cutoff) == Symbol && !(low_k_cutoff in [:auto, :none])
        throw(ArgumentError("low_k_cutoff must be a numerical value or one of [:auto, :none]"))
    elseif typeof(low_k_cutoff) <: Number && low_k_cutoff < 0
        throw(ArgumentError("low_k_cutoff must be a positive number or one of [:auto, :none]"))
    end

    # Check that timespan is compatible with chunkstep (if chunkwise solution requested)
    if solve_chunks
        try
            n_chunks_reqd = Int(tspan[2] / solve_chunkstep)
        catch e
            if e isa InexactError
                throw(ArgumentError("Simulation timespan is not divisible by requested chunkwise simulation step size"))
            else
                rethrow(e)
            end
        end
    end

    # Ensure save interval is less than chunkstep
    if solve_chunks && !isnothing(save_interval) && save_interval > solve_chunkstep
        throw(ArgumentError("Solution save interval must be less than chunkwise simulation step size"))
    end

    return ODESimulationParams(tspan, u0, solver, jac, sparse, abstol, reltol,
        adaptive_tols, update_tols, solve_chunks, solve_chunkstep,
        maxiters, ban_negatives, progress, save_interval,
        low_k_cutoff, low_k_maxconc, allow_short_u0)
end