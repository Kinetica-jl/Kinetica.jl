"""
    max_rates = get_max_rates(conditions, calculator)

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
    initial_rates = get_initial_rates(conditions, calculator)

Calculate initial rate constants for a variable condition simulation.
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
    k = calculate_discrete_rates(conditions, calculator, nr[, uType])

Calculates rate constants over a set of variable conditions.

Returns a 
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

    return ODESolution{uType, 2}(
        k_precalc,
        nothing,
        nothing,
        tstops,
        nothing,
        DummyODEProblem(; u0=k_precalc[1], tspan=[tstops[begin], tstops[end]], syms=[Symbol("k$i") for i in 1:nr]),
        nothing,
        SciMLBase.LinearInterpolation(tstops, k_precalc),
        false,
        0,
        nothing,
        nothing,
        ReturnCode.Default
    )
end


"""
    insert_inert!(rd, sd, inert_species)

Inserts inert species into all unimolecular reactions.

Updates `sd` to have fragment data for the new species,
then converts all unimolecular reactions to 'bimolecular' reactions
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
            xyz = xyz_to_frames(pbmol.write("xyz"))[1]

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
    uni_reactions = [i for i in 1:rd.nr if length(rd.reacs[i]) == 1 && rd.stoic_reacs[i][1] == 1]

    # Modify all unimolecular reactions to be bimolecular with the inert species as a collision partner.
    for (i, (species, sid)) in enumerate(zip(inert_species, inert_species_ids))
        # On all inert species except the last, create new bimolecular reactions.
        # Leave the original unimolecular reaction as a template to copy from.
        if i < length(inert_species)
            for rid in uni_reactions
                new_reacs = vcat(rd.reacs[rid], species)
                new_prods = vcat(rd.prods[rid], species)
                new_id_reacs = vcat(rd.id_reacs[rid], sid)
                new_id_prods = vcat(rd.id_prods[rid], sid)
                new_stoic_reacs = vcat(rd.stoic_reacs[rid], 1)
                new_stoic_prods = vcat(rd.stoic_prods[rid], 1)
                new_dH = rd.dH[rid]

                all_reacs = sort(reduce(vcat, [[spec for _ in new_stoic_reacs[spos]] for (spos, spec) in enumerate(new_reacs)]))
                all_prods = sort(reduce(vcat, [[spec for _ in new_stoic_prods[spos]] for (spos, spec) in enumerate(new_prods)]))
                new_rhash = stable_hash(vcat(all_reacs, all_prods))

                push!(rd.reacs, new_reacs)
                push!(rd.prods, new_prods)
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
                push!(rd.reacs[rid], species)
                push!(rd.prods[rid], species)
                push!(rd.id_reacs[rid], sid)
                push!(rd.id_prods[rid], sid)
                push!(rd.stoic_reacs[rid], 1)
                push!(rd.stoic_prods[rid], 1)

                # Update reaction hash for new reactants.
                all_reacs = sort(reduce(vcat, [[spec for _ in rd.stoic_reacs[rid][spos]] for (spos, spec) in enumerate(rd.reacs[rid])]))
                all_prods = sort(reduce(vcat, [[spec for _ in rd.stoic_prods[rid][spos]] for (spos, spec) in enumerate(rd.prods[rid])]))
                rhash = stable_hash(vcat(all_reacs, all_prods))
                rd.rhash[rid] = rhash
            end
        end
    end
end


"""
    n_removed = apply_low_k_cutoff!(rd, calc, pars, conditions)

Removes low-rate reactions from `rd` and `calc` according to the cutoff in `pars.low_k_cutoff`.

If the cutoff is a numeric value, it is used directly. If it is
`:auto`, automatically decides on a safe value where the removed
reactions would not contribute to the network over the timespan
of the simulation. If it is `:none`, does not apply a cutoff and
returns.
"""
function apply_low_k_cutoff!(rd::RxData{iType, fType}, calc::cType, 
        pars::ODESimulationParams, conditions::ConditionSet) where {
        iType, fType, cType <: AbstractKineticCalculator}

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

    # Calculate maximum rate constants.
    max_rates = get_max_rates(conditions, calc)
    
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
    make_rs(k, spec, t, rd[, combinatoric_ratelaws])

Makes a Catalyst ReactionSystem from all currently implemented reactions.

Should always be preceded by a call to `@parameters` to refresh `k` and
a call to `@variables` to refresh `species` and `t`.
"""
function make_rs(k, spec, t, rd::RxData; combinatoric_ratelaws=false)
    rxs = []
    for i in 1:rd.nr
        sr = rd.stoic_reacs[i]
        sp = rd.stoic_prods[i]
        rx = Reaction(k[i], [spec[rd.id_reacs[i][j]] for j in 1:length(sr)], [spec[rd.id_prods[i][j]] for j in 1:length(sp)], sr, sp)
        push!(rxs, rx)
    end

    @named rs = ReactionSystem(rxs, t, spec, k; combinatoric_ratelaws=combinatoric_ratelaws)
    
    return rs
end


"""
    adaptive_solve!(integrator, pars, solvecall_kwargs[, print_status])


"""
function adaptive_solve!(integrator, pars::ODESimulationParams, solvecall_kwargs::Dict{Symbol, Any}; print_status::Bool=false)
    if print_status 
        @info " - Solving..." 
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
                @info " - Solved!" 
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
    affect! = CompleteRateUpdateAffect(k_precalc)

Condition function for discrete rate update callback in complete timescale simulations.

For efficiency, and to avoid errors due to FP imprecision,
keeps track of the number of time stops it has made and
indexes precalculated rate constants with this, rather than
trying to index based on the current time.

NOTE: This may be a poor way of implementing this, as it is
incompatible with other callbacks in the same simulation. If
this is required, a multidimensional linear interpolator for
the rate constants may be necessary to deal with arbitrary
time stops.
"""
mutable struct CompleteRateUpdateAffect
    k_precalc::SciMLBase.AbstractODESolution
end
CompleteRateUpdateAffect(k_precalc) = return CompleteRateUpdateAffect(k_precalc, 1)
function (self::CompleteRateUpdateAffect)(integrator)
    integrator.p = self.k_precalc(integrator.t)
end


"""
    condition = ChunkwiseRateUpdateCondition(tstops_local)

Condition function for discrete rate update callback in chunkwise simulations.

Fires when local solver time is an element of `tstops_local`,
an array of timestops also on this local timescale.

In cases where ``τ_{update} > t_{loop}`` this may be empty, in 
which case the rate is not updated within this local loop.
"""
mutable struct ChunkwiseRateUpdateCondition
    tstops_local::Vector{Float64}
end
function (self::ChunkwiseRateUpdateCondition)(u, t, integrator)
    t ∈ self.tstops_local
end


"""
    affect! = ChunkwiseRateUpdateAffect(t_loop, n_loops)

Affect! function for discrete rate update callback in chunkwise simulations.

Fires on request of a `ChunkwiseRateUpdateCondition`. Calculates 
global simulation time using `t_chunk` and `n_chunks`, then 
updates rate constants using conditions interpolated from their
solutions on the global timescale.
"""
mutable struct ChunkwiseRateUpdateAffect{tType}
    t_chunk::tType
    n_chunks::Int
    k_precalc::AbstractDiffEqArray
end
function (self::ChunkwiseRateUpdateAffect)(integrator)
    t = integrator.t + (self.n_chunks*self.t_chunk)
    integrator.p = self.k_precalc(t)
end