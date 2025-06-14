"""
    get_max_rates(conditions::ConditionSet, calculator<:AbstractKineticCalculator)

Calculates the maximum rate constants for reactions under variable conditions.

Since rate expressions with arbitrary parameters can yield maximum
rate constants either when individual conditions are at their maxima,
at their minima, or a combination of both, all permutations of
minimum and maximum conditions must be explored to find the maximum
rate constants for a given simulation.

After determining which conditions are variable, constructs all
possible permutations of minimum/maximum variable conditions by
means of binary enumeration. Calculates all permutations of rate
constants by running condition permutations through `calculator`,
then selects the permutation with the greatest average rate as
the maximum rates.
"""
function get_max_rates(conditions::ConditionSet, calculator::AbstractKineticCalculator)
    static_condition_map = []
    minmax_variable_condition_map = []
    n_conditions = length(conditions.symbols)
    for i in 1:n_conditions
        if isstatic(conditions.profiles[i])
            push!(static_condition_map, Pair(conditions.symbols[i], conditions.profiles[i].value))
        else
            push!(minmax_variable_condition_map, [
                Pair(conditions.symbols[i], minimum(conditions.profiles[i])),
                Pair(conditions.symbols[i], maximum(conditions.profiles[i]))
            ])
        end
    end
    n_variable_conditions = length(minmax_variable_condition_map)
    if n_variable_conditions == 0
        max_rates = calculator(; static_condition_map...)
    else
        max_rate_condition_sets = []
        max_rate_condition_permutations = lpad.(string.(0:(2^n_variable_conditions)-1, base=2), n_variable_conditions, '0')
        for perm in max_rate_condition_permutations
            mrcs = []
            for (cond, sd) in enumerate(perm)
                d = parse(Int, sd)+1
                push!(mrcs, minmax_variable_condition_map[cond][d])
            end
            mrcs = vcat(mrcs, static_condition_map)
            push!(max_rate_condition_sets, mrcs)
        end

        max_rate_permutations = [calculator(; cset...) for cset in max_rate_condition_sets]
        max_rates = max_rate_permutations[findmax([mean(perm) for perm in max_rate_permutations])[2]]
    end

    return max_rates
end


"""
    get_initial_rates(conditions::ConditionSet, calculator<:AbstractKineticCalculator)

Calculates initial rate constants for a variable condition simulation.
"""
function get_initial_rates(conditions::ConditionSet, calculator::AbstractKineticCalculator)
    bound_conditions = []
    for (sym, profile) in zip(conditions.symbols, conditions.profiles)
        if isstatic(profile)
            push!(bound_conditions, Pair(sym, profile.value))
        else
            push!(bound_conditions, Pair(sym, profile.X_start))
        end
    end
    initial_rates = calculator(; bound_conditions...)
    return initial_rates    
end


"""
    calculate_discrete_rates(conditions::ConditionSet, calculator<:AbstractKineticCalculator, nr::Int[, uType=Float64])

Calculates rate constants over a set of variable conditions.

`conditions` must have a valid `tstops` array to iterate over,
usually created by `create_discrete_tstops!`. For every time in
this array, calculates all rate constants based on interpolations
of variable condition profiles.

Returns a `DiffEqArray` of rate constant vectors, one for each
time point. The numeric type of this array defaults to `Float64`,
but can be modified to match species concentrations with the 
`uType` argument.
"""
function calculate_discrete_rates(conditions::ConditionSet, calculator::AbstractKineticCalculator, nr::Int; uType=Float64)
    if !conditions.discrete_updates
        throw(ErrorException("Cannot calculate discrete rates for a continuous ConditionSet."))
    end

    tstops = get_tstops(conditions)
    scs = get_static_conditions(conditions)
    vcs = get_variable_conditions(conditions)
    k_precalc = [zeros(uType, nr) for _ in tstops]
    for (i, tstop) in enumerate(tstops)
        bound_conditions = vcat(
            scs,
            [Pair(vpair.first, vpair.second(tstop)[1]) for vpair in vcs]
        )
        k_precalc[i] = calculator(; bound_conditions...)
    end

    return DiffEqArray(k_precalc, tstops)
end


"""
    insert_inert!(rd::RxData, sd::SpeciesData, inert_species::Vector{String})

Inserts inert species into all unimolecular reactions.

Updates `sd` to have fragment data for the new species,
then converts all unimolecular reactions to bimolecular reactions
where the inert species is a bystander for the sake of being a
collision partner.

When multiple `inert_species` are present, creates new reactions
to allow for multiple reactive channels through different collision
partners.
"""
function insert_inert!(rd::RxData, sd::SpeciesData, inert_species::Vector{String})
    # Add inert species to SpeciesData, if not already present.
    inert_species_ids = []
    for species in inert_species
        if !(species in keys(sd.toInt))
            pbmol = pybel.readstring("smi", species)
            pbmol.addh()
            pbmol.make3D()
            xyz = xyz_to_frame(pbmol.write("xyz"))

            inert_id = sd.n + 1
            push!(inert_species_ids, inert_id)
            sd.toInt[species] = inert_id
            sd.toStr[inert_id] = species
            sd.xyz[inert_id] = xyz
            sd.n += 1
        else
            push!(inert_species_ids, sd.toInt[species])
        end
    end

    # Identify all unimolecular reactions.
    uni_reactions = [i for i in 1:rd.nr if length(rd.id_reacs[i]) == 1 && rd.stoic_reacs[i][1] == 1]

    # Modify all unimolecular reactions to be bimolecular with the inert species as a collision partner.
    for (i, (species, sid)) in enumerate(zip(inert_species, inert_species_ids))
        # On all inert species except the last, create new bimolecular reactions.
        # Leave the original unimolecular reaction as a template to copy from.
        if i < length(inert_species)
            for rid in uni_reactions
                new_id_reacs = vcat(rd.id_reacs[rid], sid)
                new_id_prods = vcat(rd.id_prods[rid], sid)
                new_stoic_reacs = vcat(rd.stoic_reacs[rid], 1)
                new_stoic_prods = vcat(rd.stoic_prods[rid], 1)
                new_dH = rd.dH[rid]

                all_reacs = sort(reduce(vcat, [[sd.toStr[sid] for _ in 1:new_stoic_reacs[spos]] for (spos, sid) in enumerate(rd.id_reacs[rid])]))
                all_prods = sort(reduce(vcat, [[sd.toStr[sid] for _ in 1:new_stoic_prods[spos]] for (spos, sid) in enumerate(rd.id_prods[rid])]))
                new_rhash = stable_hash(vcat(all_reacs, all_prods); version=4)

                # Doesn't make sense to update the atom maps with inert species since they don't really exist...
                push!(rd.mapped_rxns, rd.mapped_rxns[rid])
                push!(rd.id_reacs, new_id_reacs)
                push!(rd.id_prods, new_id_prods)
                push!(rd.stoic_reacs, new_stoic_reacs)
                push!(rd.stoic_prods, new_stoic_prods)
                push!(rd.dH, new_dH)
                push!(rd.rhash, new_rhash)
                rd.nr += 1
            end
        # On the final inert species, modify the unimolecular template reaction.
        else
            for rid in uni_reactions
                push!(rd.id_reacs[rid], sid)
                push!(rd.id_prods[rid], sid)
                push!(rd.stoic_reacs[rid], 1)
                push!(rd.stoic_prods[rid], 1)

                # Update reaction hash for new reactants.
                all_reacs = sort(reduce(vcat, [[sd.toStr[sid] for _ in 1:rd.stoic_reacs[rid][spos]] for (spos, sid) in enumerate(rd.id_reacs[rid])]))
                all_prods = sort(reduce(vcat, [[sd.toStr[sid] for _ in 1:rd.stoic_prods[rid][spos]] for (spos, sid) in enumerate(rd.id_prods[rid])]))
                rhash = stable_hash(vcat(all_reacs, all_prods); version=4)
                rd.rhash[rid] = rhash
            end
        end
    end
end


"""
    apply_low_k_cutoff!(rd::RxData, calc<:AbstractKineticCalculator, pars::ODESimulationParams, conditions::ConditionSet)

Removes low-rate reactions from `rd` and `calc` according to the cutoff in `pars.low_k_cutoff`.

If the cutoff is a numeric value, it is used directly. If it is
`:auto`, automatically decides on a safe value where the removed
reactions would not contribute to the network over the timespan
of the simulation. If it is `:none`, does not apply a cutoff.

Multiplies the calculated maximum rates of reaction by a theoretical
maximum concentration that any reactants could attain through a 
kinetic simulation, as set by `pars.low_k_maxconc`. This concentration
multiplier is squared to emulate a bimolecular reaction.

Returns the number of low-rate reactions that have been removed
from the CRN.
"""
function apply_low_k_cutoff!(rd::RxData{iType, fType}, calc::cType, 
        pars::ODESimulationParams, conditions::ConditionSet) where {
        iType, fType, cType <: AbstractKineticCalculator}

    # Establish what the value of the cutoff should be.
    if pars.low_k_cutoff == :none
        @info "   - Low rate cutoff: none"
        return 0
    elseif pars.low_k_cutoff == :auto
        k_cutoff = pars.reltol/pars.tspan[end]
        @info "   - Low rate cutoff: automatic (cutoff = $(k_cutoff))"
    else
        k_cutoff = uType(pars.low_k_cutoff)
        @info "   - Low rate cutoff: manual (cutoff = $(k_cutoff))"
    end

    # Calculate maximum rate constants.
    max_rates = get_max_rates(conditions, calc) .* pars.low_k_maxconc^2
    
    # Find low rate reactions.
    low_rate_ids = Vector{iType}()
    for (i, rate) in enumerate(max_rates)
        if rate < k_cutoff
            push!(low_rate_ids, i)
        end
    end

    # Remove from network.
    splice!(rd, calc, low_rate_ids)

    @info "   - Removed $(length(low_rate_ids)) low-rate reactions from network."
    return length(low_rate_ids)
end


"""
    make_u0(sd::SpeciesData, pars::ODESimulationParams)

Construct the initial concentration vector `u0` from the input in `pars.u0`.

Converts the input initial concentrations into a vector of length
`sd.n`. If `pars.u0` is a vector, this is by default only 
allowed if the input vector has an entry for every species. This 
behaviour can be changed with `pars.allow_short_u0`, which
fills any remaining species concentrations with zeros.

If the input initial concentrations are in a Dict, converts species
names to IDs and correctly populates an array at the right indeces.
"""
function make_u0(sd::SpeciesData, pars::ODESimulationParams)
    if pars.u0 isa Vector
        @info "   - Starting with initial concentrations for all species.."
        if length(pars.u0) != sd.n
            if pars.allow_short_u0
                @info "   - pars.u0 shorter than expected, setting all trailing species to 0"
                u0 = zeros(eltype(pars.u0), sd.n)
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
            if spec in keys(sd.toInt)
                species_ids[i] = sd.toInt[spec]
            else
                throw(ErrorException("Species $spec not in SpeciesData. Check pars.u0 is correct."))
            end
        end

        # Create populated u0 vector.
        utype = valtype(pars.u0)
        u0 = zeros(utype, sd.n)
        for id in species_ids
            u0[id] = pars.u0[sd.toStr[id]]
        end
    end
    return u0
end


"""
    make_rs(k, spec, t, rd::RxData[, variable_k=false, combinatoric_ratelaws=false])
    make_rs(k, spec, t, rd::RxData, p[, combinatoric_ratelaws=false])
    
Makes a Catalyst `ReactionSystem`` from all currently implemented reactions.

`k` should be a vector of rate constants, defined either as MTK
variables or parameters. If performing a chunkwise simulation with
`k` as variables, chunk parameters should be passed via `p`. If
performing a complete timescale simulation, `variable_k` should be
true.

`spec` should be a vector of species defined as Catalyst species.
`t` should be the MTK variable for simulation time.

If a SSA-based simulation is required, `combinatoric_ratelaws` can
be enabled. Otherwise, this should be left disabled.
"""
function make_rs(k, spec, t, rd::RxData; variable_k=false, combinatoric_ratelaws=false)
    rxs = []
    for i in 1:rd.nr
        sr = rd.stoic_reacs[i]
        sp = rd.stoic_prods[i]
        rx = Reaction(k[i], [spec[rd.id_reacs[i][j]] for j in 1:length(sr)], [spec[rd.id_prods[i][j]] for j in 1:length(sp)], sr, sp)
        push!(rxs, rx)
    end

    if variable_k
        @named rs = ReactionSystem(rxs, t; combinatoric_ratelaws=combinatoric_ratelaws)
    else
        @named rs = ReactionSystem(rxs, t, spec, k; combinatoric_ratelaws=combinatoric_ratelaws)
    end
    
    return rs
end

function make_rs(k, spec, t, rd::RxData, p; combinatoric_ratelaws=false)
    rxs = []
    for i in 1:rd.nr
        sr = rd.stoic_reacs[i]
        sp = rd.stoic_prods[i]
        rx = Reaction(k[i], [spec[rd.id_reacs[i][j]] for j in 1:length(sr)], [spec[rd.id_prods[i][j]] for j in 1:length(sp)], sr, sp)
        push!(rxs, rx)
    end

    states = vcat(spec, k)
    @named rs = ReactionSystem(rxs, t, states, p; combinatoric_ratelaws=combinatoric_ratelaws)
    
    return rs
end


"""
    adaptive_solve!(integrator, pars::ODESimulationParams, solvecall_kwargs::Dict{Symbol, Any}[, print_status=false])

Tries to solve `integrator` with increasing tolerances as solution becomes unstable.

Repeatedly attempts to call `solve!(integrator)`, starting with
the `abstol` and `reltol` defined in `pars`. If this fails,
reduces both values by an order of magnitude and tries again.

This can continue until either the tolerances cannot be reduced
without being less than the precision of the numeric type being
used within the solution, or until 5 attempts have elapsed, at
which point an exception will be thrown.

This behaviour can be disabled with `pars.adaptive_tols=false`,
which will cause this function to fail after the first failed
attempt. 

If the integrator finishes successfully and `pars.update_tols=true`,
the final solver tolerances used in this successful simulation are
written back to `pars`. These tolerances are taken from their
entries in `solvecall_kwargs`, which are always updated during
successive solution attempts.
"""
function adaptive_solve!(integrator, pars::ODESimulationParams, solvecall_kwargs::Dict{Symbol, Any}; print_status=false)
    if print_status 
        @info " - Solving network..." 
        flush_log()
    end

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
            if print_status 
                @info " - Solved!\n" 
                flush_log()
            end
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
                flush_log()
                integrator.opts.abstol = solvecall_kwargs[:abstol]
                integrator.opts.reltol = solvecall_kwargs[:reltol]
                reinit!(integrator)
            end
        end
    end
end


"""
    CompleteRateUpdateAffect(k_precalc<:AbstractDiffEqArray)

Affect! function for discrete rate update callback in complete timescale simulations.

Calculates new rate constants from direct interpolation on
`k_precalc::AbstractDiffEqArray`.
"""
mutable struct CompleteRateUpdateAffect{idxType}
    k_precalc::SciMLBase.AbstractDiffEqArray
    k_idx_map::Vector{idxType}
end

"""
    (self::CompleteRateUpdateAffect)(integrator)

Sets rate constants at `integrator.p` to those in `self.k_precalc` at time `integrator.t`.
"""
function (self::CompleteRateUpdateAffect)(integrator)
    k_vals = self.k_precalc(integrator.t)
    for i in 1:length(self.k_idx_map)
        setindex!(integrator.p, k_vals[i], self.k_idx_map[i])
    end
end


"""
    ChunkwiseRateUpdateCondition(tstops_local::Vector{Float64})

Condition function for discrete rate update callback in chunkwise simulations.

Fires when local solver time is an element of `tstops_local`,
an array of timestops also on this local timescale.

In cases where ``τ_{update} > t_{loop}`` this may be empty, in 
which case the rate is not updated within this local loop.
"""
mutable struct ChunkwiseRateUpdateCondition
    tstops_local::Vector{Float64}
end

"""
    (self::ChunkwiseRateUpdateCondition)(u, t, integrator)

Checks if current time `t` is within `self.tstops_local`.
"""
function (self::ChunkwiseRateUpdateCondition)(u, t, integrator)
    t ∈ self.tstops_local
end


"""
    ChunkwiseRateUpdateAffect(t_chunk, n_chunks::Int, k_precalc<:AbstractDiffEqArray)

Affect! function for discrete rate update callback in chunkwise simulations.

Fires on request of a `ChunkwiseRateUpdateCondition`. Calculates 
global simulation time using `t_chunk` and `n_chunks`, then 
updates rate constants using conditions interpolated from their
solutions on the global timescale.
"""
mutable struct ChunkwiseRateUpdateAffect{tType, idxType}
    t_chunk::tType
    n_chunks::Int
    k_precalc::AbstractDiffEqArray
    k_idx_map::Vector{idxType}
end

"""
    (self::ChunkwiseRateUpdateAffect)(integrator)

Sets rate constants at `integrator.p` to those in `self.k_precalc` at current local chunk time.

Global simulation time is calculated from addition of local 
chunk time and number of previous chunks.
"""
function (self::ChunkwiseRateUpdateAffect)(integrator)
    t = integrator.t + (self.n_chunks*self.t_chunk)
    k_vals = self.k_precalc(t)
    for i in 1:length(self.k_idx_map)
        setindex!(integrator.p, k_vals[i], self.k_idx_map[i])
    end
end