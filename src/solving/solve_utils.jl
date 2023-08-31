"""
    n_removed = apply_low_k_cutoff!(rd, pars, rates, cutoff)

Removes low-rate reactions from `rd` according to the cutoff in `pars.low_k_cutoff`.

If the cutoff is a numeric value, it is used directly. If it is
`:auto`, automatically decides on a safe value where the removed
reactions would not contribute to the network over the timespan
of the simulation. If it is `:none`, does not apply a cutoff and
returns.
"""
function apply_low_k_cutoff!(rd::RxData{iType, fType}, pars::ODESimulationParams,
        rates::Vector{uType}) where {iType, fType, uType <: Base.AbstractFloat}
    # Establish what the value of the cutoff should be.
    if pars.low_k_cutoff == :none
        @info " - Low rate cutoff: none"
        return 0
    elseif pars.low_k_cutoff == :auto
        k_cutoff = pars.reltol/pars.tspan[end]
        @info " - Low rate cutoff: automatic (cutoff = $(k_cutoff))"
    else
        k_cutoff = uType(pars.low_k_cutoff)
        @info " - Low rate cutoff: manual (cutoff = $(k_cutoff))"
    end
    
    # Find low rate reactions.
    low_rate_ids = Vector{iType}()
    for (i, rate) in enumerate(rates)
        if rate < k_cutoff
            push!(low_rate_ids, i)
        end
    end

    # Remove from network.
    if length(low_rate_ids) > 0
        for f in fieldnames(typeof(rd))
            if f != :nr
                splice!(getfield(rd, f), low_rate_ids)
            end
        end
        rd.nr -= length(low_rate_ids)
        splice!(rates, low_rate_ids)
    end

    @info "   - Removed $(length(low_rate_ids)) low-rate reactions from network."
    return length(low_rate_ids)
end


"""
    u0 = make_u0(species, pars)

Construct the initial concentration vector `u0` from the input in `pars.u0`.

Converts the input initial concentrations into a vector of length
`species.n`. If `pars.u0` is a vector, this is by default only 
allowed if the input vector has an entry for every species. This 
behaviour can be changed with `pars.allow_short_u0`, which
fills any remaining species concentrations with zeros.

If the input initial concentrations are in a Dict, converts species
names to IDs and correctly populates an array at the right indeces.
"""
function make_u0(species::SpeciesData, pars::ODESimulationParams)
    if pars.u0 isa Vector
        @info "   - Starting with initial concentrations for all species.."
        if length(pars.u0) != species.n
            if pars.allow_short_u0
                @info "   - pars.u0 shorter than expected, setting all trailing species to 0"
                u0 = zeros(eltype(pars.u0), species.n)
                u0[1:length(pars.u0)] = pars.u0
            else
                throw(ErrorException("Length of supplied initial concentration vector does not match with number of species in system."))
            end
        else
            u0 = pars.u0
        end
    # If a dict, ensure it is pointing to the correct species.
    else
        initial_species = keys(pars.u0)
        # Check requested species exist in known fragments.
        species_ids = zeros(Int, length(initial_species))
        for (i, spec) in enumerate(initial_species)
            if spec in keys(species.toInt)
                species_ids[i] = species.toInt[spec]
            else
                throw(ErrorException("Species $spec not in SpeciesData. Check pars.u0 is correct."))
            end
        end

        # Create populated u0 vector.
        utype = valtype(pars.u0)
        u0 = zeros(utype, species.n)
        for id in species_ids
            u0[id] = pars.u0[species.toStr[id]]
        end
    end
    return u0
end


"""
    make_rs(k, spec, t, rd[, constraints, combinatoric_ratelaws])

Makes a Catalyst ReactionSystem from all currently implemented reactions.

Should always be preceded by a call to `@parameters` to refresh `k` and
a call to `@variables` to refresh `species` and `t`.
"""
function make_rs(k, spec, t, rd::RxData; constraints=nothing, combinatoric_ratelaws=false)
    rxs = []
    for i in 1:rd.nr
        sr = rd.stoic_reacs[i]
        sp = rd.stoic_prods[i]
        rx = Reaction(k[i], [spec[rd.id_reacs[i][j]] for j in 1:length(sr)], [spec[rd.id_prods[i][j]] for j in 1:length(sp)], sr, sp)
        push!(rxs, rx)
    end

    if isnothing(constraints)
        @named rs = ReactionSystem(rxs, t, collect(spec), collect(k); 
                                   combinatoric_ratelaws=combinatoric_ratelaws)
    else
        @named rs = ReactionSystem(rxs, t, collect(spec), collect(k); 
                                   combinatoric_ratelaws=combinatoric_ratelaws, 
                                   checks=false, systems=[constraints])
    end
    
    return rs
end


"""
    adaptive_solve!(integrator, pars, solvecall_kwargs[, print_status])


"""
function adaptive_solve!(integrator, pars::ODESimulationParams, solvecall_kwargs::Dict{Symbol, Any}; print_status::Bool=false)
    if print_status @info " - Solving..." end

    # Attempt solutions with lower tolerances to attempt to ensure convergence.
    success = false
    iters = 0
    mintol = eps(eltype(integrator.sol.prob.u0))
    while !success
        iters += 1
        with_global_logger() do
            solve!(integrator)
        end
        if SciMLBase.successful_retcode(integrator.sol.retcode)
            success = true
            if print_status @info " - Solved!" end
            if pars.update_tols && solvecall_kwargs[:abstol] != opars.abstol
                @info "   - Writing new tolerances back to ODEParams."
                pars.abstol = solvecall_kwargs[:abstol]
                pars.reltol = solvecall_kwargs[:reltol]
            end
        else
            if !pars.adaptive_tols
                @error " - Solve failed, not retrying as adaptive tolerance is not enabled."
                throw(ErrorException("ODE solution failed."))
            elseif iters >= 5
                @error " - Too many attempts have been made to reduce solver tolerance, exiting."
                throw(ErrorException("ODE solution failed."))
            elseif solvecall_kwargs[:abstol]/10 <= mintol || solvecall_kwargs[:reltol]/10 <= mintol
                @error " - Solution cannot be converged by reducing solver tolerance any further, exiting."
                throw(ErrorException("ODE solution failed."))
            else
                solvecall_kwargs[:abstol] /= 10
                solvecall_kwargs[:reltol] /= 10
                @warn "   - ODE solution failed at current solver tolerances."
                @warn "   - Reducing tolerances to abstol = $(solvecall_kwargs[:abstol]) reltol = $(solvecall_kwargs[:reltol])"
                integrator.opts.abstol = solvecall_kwargs[:abstol]
                integrator.opts.reltol = solvecall_kwargs[:reltol]
                reinit!(integrator)
            end
        end
    end
end