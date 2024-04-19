abstract type AbstractExploreMethod end

"""
    DirectExplore(rdir_head::String, reac_smiles::VEctor{String}, cde::CDE[, maxiters::Int=1000, rxn_convergence_threshold::Int=5])

Keyword-based container for parameters used in direct CRN exploration.

Contains fields for:
* Top level of CRN exploration directory (`rdir_head`)
* SMILES string(s) of main breakdown reactant(s) being studied (`reac_smiles`)
* `CDE` instance (`cde`)
* Maximum number of iterations to perform (`maxiters`)
* Number of iterations with no change in reactions to consider as converged (`rxn_convergence_threshold`)
"""
@kwdef mutable struct DirectExplore <: AbstractExploreMethod
    rdir_head::String
    reac_smiles::Vector{String}
    cde::CDE
    maxiters::Integer = 1000
    rxn_convergence_threshold::Integer = 5
end


"""
    IterativeExplore(rdir_head::String, reac_smiles::VEctor{String}, cde::CDE[, maxiters::Int=1000, 
                     rxn_convergence_threshold::Int=5, seed_convergence_threshold::Int=3, seed_conc=0.05,
                     independent_blacklist::Vector{String}=[], inert_species::Vector{String}=[]])

Keyword-based container for parameters used in iterative kinetics-based CRN exploration.

Contains fields for:
* Top level of CRN exploration directory (`rdir_head`)
* SMILES string(s) of main breakdown reactant(s) being studied (`reac_smiles`)
* `CDE` instance (`cde`)
* Maximum number of iterations to perform (`maxiters`)
* Number of subspace iterations with no change in reactions to consider a subspace converged (`rxn_convergence_threshold`)
* Number of level iterations with no change in seeds to consider the network converged (`seed_convergence_threshold`)
* Concentration above which species will be selected as seeds each level (`seed_conc`)
* Blacklist of species to avoid doing independent subspace explorations on (`independent_blacklist`)
* Inert species that should not be considered for reaction (`inert_species`)
"""
@kwdef struct IterativeExplore{uType} <: AbstractExploreMethod
    rdir_head::String
    reac_smiles::Vector{String}
    cde::CDE
    maxiters::Integer = 1000
    rxn_convergence_threshold::Integer = 5
    seed_convergence_threshold::Integer = 3
    seed_conc::uType = 0.05
    independent_blacklist::Vector{String} = String[]
    inert_species::Vector{String} = String[]
end


"""
    explore_network(exploremethod::DirectExplore, solvemethod[, savedir])
    explore_network(exploremethod::IterativeExplore, solvemethod[, savedir])

Runs network exploration with one of the available methods.

If `exploremethod isa DirectExplore`, runs a single-level
network exploration to attempt to locate all relevant
reactions in a radius of `exploremethod.cde.radius`
species from the starting system.

If `exploremethod isa IterativeExplore`, runs a multi-level
iterative network exploration using kinetic simulations to
identify seed species for successive levels, in order to
fully characterise the reaction space relevant to the
conditions in `solvemethod.conditions`.
"""
function explore_network end

function explore_network(exploremethod::DirectExplore,
                         solvemethod::AbstractSolveMethod;
                         savedir::Union{String, Nothing}=nothing)

    @info "#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#"
    @info "Kinetica Direct CRN Exploration"
    @info "#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#\n"
    flush_log()

    if !isdir(exploremethod.rdir_head) mkdir(exploremethod.rdir_head) end
    if !isnothing(savedir) && !isdir(savedir) mkdir(savedir) end

    loc = find_current_loc(exploremethod.rdir_head)
    if loc.level > 1
        error("Current CRN level is greater than 1. Are you trying to continue an iterative exploration?")
    end
    seeds = exploremethod.reac_smiles
    if loc.level == 0
        sd, rd = init_network()
        for rsmi in unique(seeds)
            # Ensure initial reactants are consistently across CRN generations when RNG is seeded.
            rxyz = xyz_to_frame(xyz_from_smiles(rsmi, generator=:rdkit, seed=rand(1:999999999)))
            push_unique!(sd, rsmi, rxyz)
        end
        inc_level!(loc)
        setup_level(loc, sd, seeds)
        @info "Starting breakdown generation within a radius of $(exploremethod.cde.radius) reactions.\n"
    else
        cleanup_network(loc.rdir_head)
        sd, rd = import_network(loc.rdir_head)
        @info "Continuing breakdown generation within a radius of $(exploremethod.cde.radius) reactions.\n"
    end
    
    n_seeds = length(seeds)
    n_subspaces = n_seeds == 1 ? 1 : n_seeds + 1
    explored_seeds = []
    while loc.subspace < n_subspaces
        spec = seeds[loc.subspace]
        if spec in explored_seeds
            cpath = joinpath(pathof(loc), "isconv")
            open(cpath, "w") do f
                write(f, "true")
            end
            @info "Same-species reactions between $(spec) already covered."
            @info "Skipping Subspace $(loc.subspace)\n"; flush_log()
        else
            explore_subspace!(sd, rd, loc, exploremethod)
            push!(explored_seeds, spec)
        end
        inc_subspace!(loc)
    end

    explore_subspace!(sd, rd, loc, exploremethod)
    @info "Exploration complete, running kinetic simulation of current network."; flush_log()
    res = solve_network(solvemethod, sd, rd)
    @info "Direct network exploration complete."

    if !isnothing(savedir)
        @info "Saving finished network."
        saveto = joinpath(savedir, "direct_network_final.bson")
        save_output(res, saveto)
        @info "Network saved to $saveto"
    end

    return res
end

function explore_network(exploremethod::IterativeExplore,
                         solvemethod::AbstractSolveMethod;
                         savedir::Union{String, Nothing}=nothing)

    @info "#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-"
    @info "Kinetica Iterative CRN Exploration"
    @info "#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-\n"
    flush_log()

    if !isdir(exploremethod.rdir_head) mkdir(exploremethod.rdir_head) end
    if !isnothing(savedir) && !isdir(savedir) mkdir(savedir) end

    # Determine if continuing an existing exploration or starting from scratch.
    loc = find_current_loc(exploremethod.rdir_head); flush_log()
    if loc.level == 0
        sd, rd = init_network()
        make_inert_file(exploremethod.rdir_head, exploremethod.inert_species)
        for rsmi in vcat(exploremethod.reac_smiles, exploremethod.inert_species)
            rxyz = xyz_to_frame(xyz_from_smiles(rsmi, generator=:rdkit, seed=rand(1:999999999)))
            push_unique!(sd, rsmi, rxyz)
        end
        explored_seeds = String[]
        current_seeds = exploremethod.reac_smiles
        inc_level!(loc)
    else
        cleanup_network(loc.rdir_head)
        sd, rd = import_network(loc.rdir_head)
        explored_seeds = load_past_seeds(loc)
        current_seeds = load_current_seeds(loc)
    end

    # Loop over levels of reaction space until converged
    do_next_level = true
    convergence_count = 0
    res = nothing
    while do_next_level
        @info "##########################"
        @info "ENTERING LEVEL $(loc.level)"
        @info "##########################\n"
        flush_log()

        # Create directory structure for level (if not already created)
        setup_level(loc, sd, current_seeds)
        n_seeds = length(current_seeds)
        n_subspaces = n_seeds == 1 ? 1 : n_seeds + 1

        @info "Exploring reaction subspaces ($(loc.subspace)/$(n_subspaces)).\n"; flush_log()
        while loc.subspace < n_subspaces
            spec = current_seeds[loc.subspace]
            if spec in explored_seeds
                cpath = joinpath(pathof(loc), "isconv")
                open(cpath, "w") do f
                    write(f, "true")
                end
                @info "Same-species reactions between $(spec) already covered in prior levels."
                @info "Skipping Subspace $(loc.subspace)\n"; flush_log()
            elseif spec in exploremethod.independent_blacklist
                cpath = joinpath(pathof(loc), "isconv")
                open(cpath, "w") do f
                    write(f, "true")
                end
                @info "Same-species reactions between $(spec) prohibited by blacklist."
                @info "Skipping Subspace $(loc.subspace)\n"; flush_log()
            else
                explore_subspace!(sd, rd, loc, exploremethod)
            end
            inc_subspace!(loc)
        end

        explore_subspace!(sd, rd, loc, exploremethod)
        @info "Exploration complete, running kinetic simulation of current network."; flush_log()
        res = solve_network(solvemethod, sd, rd)

        if !isnothing(savedir)
            @info "Saving incomplete network..."; flush_log()
            saveto = joinpath(savedir, "level_network_1-$(loc.level).bson")
            save_output(res, saveto)
            @info "Network saved to $saveto\n"; flush_log()
        end

        for seed in current_seeds push!(explored_seeds, seed) end
        seeds_out_saveto = isnothing(savedir) ? nothing : joinpath(savedir, "seeds_level$(loc.level).out")
        next_seeds = identify_next_seeds(res.sol, sd, exploremethod.seed_conc;
                                        ignore=exploremethod.inert_species,
                                        saveto=seeds_out_saveto)

        if Set(current_seeds) == Set(next_seeds)
            convergence_count += 1
            if convergence_count >= exploremethod.seed_convergence_threshold
                @info "##########################"
                @info "NO NEW SEEDS FOUND FOR $(convergence_count)/$(exploremethod.seed_convergence_threshold) LEVELS"
                @info "ITERATIVE EXPLORATION COMPLETE"
                @info "##########################"
                do_next_level = false
            else
                @info "No new seeds found for $(convergence_count)/$(exploremethod.seed_convergence_threshold) levels."
                @info "Continuing to next level."
                inc_level!(loc)
                reset_subspace!(loc)
            end
        else
            @info "New seeds found, continuing to next level."
            inc_level!(loc)
            reset_subspace!(loc)
        end
        current_seeds = deepcopy(next_seeds)
        flush_log()
    end

    return res
end


"""
    explore_subspace!(sd::SpeciesData, rd::RxData, loc::ExploreLoc, exploremethod<:AbstractExploreMethod)

Finds reactions stemming from the seeds in the subspace defined by `loc`.

Uses CDE to generate mechanisms of randomly sampled reactions
between subspace seeds. Generation stops forcibly with an error
if `exploremethod.maxiters` iterations are reached. Generation
completes successfully when `exploremethod.rxn_convergence_threshold` 
iterations have passed without any new reactions being added to
the network. 
"""
function explore_subspace!(sd::SpeciesData, rd::RxData, loc::ExploreLoc, exploremethod::AbstractExploreMethod)
    @info "--------------------------"
    @info "ENTERING SUBSPACE $(loc.subspace)"
    @info "--------------------------\n"
    flush_log()

    # Check if subspace is already converged.
    cpath = joinpath(pathof(loc), "isconv")
    if isfile(cpath)
        @info "Subspace is already converged.\n"
        return
    end

    # Set CDE to work in this subspace's directory.
    exploremethod.cde.rdir = pathof(loc)
    exploremethod.cde.init_xyz = joinpath(pathof(loc), "seeds.xyz")

    # Establish where we are in this subspace.
    rcount = make_rcount(joinpath(pathof(loc), "rcount"))

    # Begin looping over chemical reaction space.
    isconv = false
    counter = 0
    no_new_reacs_iters = 0
    if rcount == 0
        @info " - Starting iterations.\n"
    else
        @info " - Continuing iterations.\n"
    end
    flush_log()

    while !isconv
        # Check we aren't iterating to infinity.
        if counter > exploremethod.maxiters
            error("$(exploremethod.maxiters) iterations exceeded, exiting loop")
        end

        counter += 1
        @info "--- ITERATION $counter ---"
        rcount += 1
        rcountrange = nothing

        # Run CDE to generate a new mechanism(s)
        if exploremethod.cde.parallel_runs > 1
            rcountrange = rcount:rcount+exploremethod.cde.parallel_runs-1
            rcountend = exploremethod.cde(rcountrange)
            if rcountend < rcountrange.start
                @warn "Sampling failed, cycling..."; flush_log()
                rcount -= 1
                continue
            end
            rcountrange = rcountrange.start:rcountend
        else
            cde_success = exploremethod.cde(rcount)
            if !cde_success
                @warn "Sampling failed, cycling..."; flush_log()
                rcount -= 1
                continue
            end
        end

        # Load in new mechanism data.
        @info " - Importing generated reactions."
        n_reacs_prev = rd.nr
        if exploremethod.cde.parallel_runs > 1
            for rc in rcountrange
                import_mechanism!(sd, rd, pathof(loc), rc)
            end
            rcount += length(rcountrange)-1
        else
            import_mechanism!(sd, rd, pathof(loc), rcount)
        end
        @info "   - Reaction network now contains $(rd.nr) reactions over $(sd.n) unique fragments."
        flush_log()

        # Check if rd has changed, increment convergence counter.
        # Rationale behind this is that the network will not be converged if new
        # reactions are still being discovered.
        if n_reacs_prev != rd.nr
            no_new_reacs_iters = 0
            @info " - New reactions discovered, reaction network not converged."
            @info " - Skipping to next iteration.\n"
            flush_log()
            continue
        end

        no_new_reacs_iters += 1
        @info " - No new reactions discovered for $(no_new_reacs_iters)/$(exploremethod.rxn_convergence_threshold) iterations.\n"

        # Check for convergence of reaction network.
        if no_new_reacs_iters >= exploremethod.rxn_convergence_threshold
            @info "   - Species subspace converged!\n"
            cpath = joinpath(pathof(loc), "isconv")
            open(cpath, "w") do f
                write(f, "true")
            end
            isconv = true
        else
            continue
        end
    end
    flush_log()
    return
end

