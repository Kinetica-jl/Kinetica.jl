"""
    rcount = make_rcount(path)

Determines the number of reaction mechanisms generated.

Reads the number of generated reactions from the file at
`path`. If this file does not exist, creates it and sets
`rcount` to zero.
"""
function make_rcount(path::String)
    rcount = 0
    if isfile(path)
        open(path, "r") do f
            rcount = parse(Int, readline(f))
        end
    else
        open(path, "w") do f
            write(f, "00000")
        end
    end

    return rcount
end


"""
    make_inert_file(path, inert_species)

Makes an inert species file at the given directory `dir`.

If an empty vector of species is supplied, returns without creating
a file. Otherwise, creates a file listing all inert species, overwriting
any other existing species file.
"""
function make_inert_file(dir::String, inert_species::Vector{String})
    if length(inert_species) == 0
        return
    else
        open(joinpath(dir, "inert.in"), "w") do f
            for spec in inert_species
                write(f, "$spec\n")
            end
        end
    end
    return
end


"""
    sd, rd = import_mechanism(rdir, rcount[, max_molecularity])

Create a CRN's initial `SpeciesData` and `RxData` from a CDE generated mechanism(s).
"""
function import_mechanism(rdir::String, rcount; max_molecularity=2)
    rsmis, rxyzs, psmis, pxyzs, dHs = ingest_cde_run(rdir, rcount)
    all_smis = vcat(reduce(vcat, rsmis), reduce(vcat, psmis))
    all_xyzs = vcat(reduce(vcat, rxyzs), reduce(vcat, pxyzs))
    sd = SpeciesData(all_smis, all_xyzs)
    rd = RxData(sd, rsmis, psmis, dHs; max_molecularity=max_molecularity)
    return sd, rd
end

"""
    import_mechanism!(sd, rd, rdir, rcount[, max_molecularity])

Extend a CRN's `SpeciesData` and `RxData` from a CDE generated mechanism(s).
"""
function import_mechanism!(sd::SpeciesData, rd::RxData, rdir::String, rcount;
        max_molecularity=2)
    rsmis, rxyzs, psmis, pxyzs, dHs = ingest_cde_run(rdir, rcount)
    all_smis = vcat(reduce(vcat, rsmis), reduce(vcat, psmis))
    all_xyzs = vcat(reduce(vcat, rxyzs), reduce(vcat, pxyzs))
    push_unique!(sd, all_smis, all_xyzs)
    push!(rd, sd, rsmis, psmis, dHs; max_molecularity=max_molecularity)
    return
end


"""
    sd, rd = import_network(rdir_head)

Imports a network from a level tree within `rdir_head`.

Recurses through level directories in `rdir_head`, then recurses
through subspace surectories in each level. Within each subspace,
imports all mechanisms within each generated reaction.

Returns the resulting network (an instance of `SpeciesData` and 
`RxData`).
"""
function import_network(rdir_head::String)
    @info "Importing all reactions in level tree under $(rdir_head)"; flush_log()
    level_dirs = readdir(rdir_head)
    level_dirs = level_dirs[startswith.(level_dirs, "level_")]
    if length(level_dirs) == 0
        error("ERROR: No network levels found in rdir_head.")
    end

    inert_file = joinpath(rdir_head, "inert.in")
    inert_species = []
    if isfile(inert_file)
        @debug "Inert species file found."
        inert_species = [line for line in eachline(inert_file) if length(line) > 0]
        @debug "Loaded inert species: $(inert_species)"
    end
    
    # Initialise network with blank sd and rd.
    sd, rd = init_network()

    # Add inert species, if there are any.
    for spec in inert_species
        mol = mol_from_smiles(spec)
        push_unique!(sd, spec, mol)
    end

    # Loop through each level, adding each subspace.
    for lv in level_dirs
        lv_dir = joinpath(rdir_head, lv)
        ss_dirs = readdir(lv_dir)
        ss_dirs = ss_dirs[startswith.(ss_dirs, "subspace_")]
        for ss in ss_dirs
            ss_dir = joinpath(lv_dir, ss)
            rcount = make_rcount(joinpath(ss_dir, "rcount"))
            for reac in 1:rcount
                import_mechanism!(sd, rd, ss_dir, reac)
            end
        end
    end

    @info "Finished network import."
    @info "Network contains $(sd.n) species over $(rd.nr) reactions.\n"
    flush_log()

    return sd, rd
end


"""
    cleanup_network(rdir_head)

Removes incomplete CDE runs from a network.

Recurses through level directories in `rdir_head` and their
respective subspace directories, checking for reaction
folders that are numbered above that subspace's `rcount`.
Removes these folders to avoid overwriting-based errors when
continuing reaction generation.
"""
function cleanup_network(rdir_head::String)
    level_dirs = readdir(rdir_head)
    level_dirs = level_dirs[startswith.(level_dirs, "level_")]
    if length(level_dirs) == 0
        return
    end

    @debug "Running cleanup to remove incomplete CDE runs."
    n_reacs_removed = 0
    for lv in level_dirs
        lv_dir = joinpath(rdir_head, lv)
        ss_dirs = readdir(lv_dir)
        ss_dirs = ss_dirs[startswith.(ss_dirs, "subspace_")]
        for ss in ss_dirs
            ss_dir = joinpath(lv_dir, ss)
            reac_dirs = readdir(ss_dir)
            reac_dirs = reac_dirs[startswith.(reac_dirs, "reac_")]
            n_reacs = length(reac_dirs)
            rcount = make_rcount(joinpath(ss_dir, "rcount"))
            for rxn in rcount+1:n_reacs
                rm(joinpath(ss_dir, reac_dirs[rxn]); recursive=true)
                n_reacs_removed += 1
            end
        end
    end

    @debug "$n_reacs_removed incomplete CDE runs removed."
end


"""
    setup_level(loc, sd, seeds)

Sets up a new level of iterative network exploration.

Uses the seed species in `seeds` to create a seed mapping
file `seeds.in` for the current level. Generates a `seeds.xyz`
for each subspace to be explores within the level, including
same-species and cross-species subspaces.
"""
function setup_level(loc::ExploreLoc, sd::SpeciesData, seeds::Vector{String})
    lvdir = pathof(loc; to_level=true)
    # Check existence of already set up level
    if isdir(lvdir)
        if isfile(joinpath(lvdir, "seeds.in"))
            @info "Level has been previously set up."
            return
        end
    else
        mkdir(lvdir)
    end

    @info "Setting up level directory tree in $lvdir"

    # Create seed mapping file
    open(joinpath(lvdir, "seeds.in"), "w") do f
        write(f, "$(length(seeds))\n")
        write(f, "SID   SMILES\n")
        for (sid, smi) in enumerate(seeds)
            write(f, "$sid    $smi\n")
        end
    end

    # Generate subspace systems for same-species reactions
    for (i, smi) in enumerate(seeds)
        ssdir = joinpath(lvdir, "subspace_$(lpad(i, 3, "0"))")
        mkdir(ssdir)
        xyz = sd.xyz[sd.toInt[smi]]
        system_from_mols([deepcopy(xyz), deepcopy(xyz)], joinpath(ssdir, "seeds.xyz"))
    end

    # Generate system for all cross-species reactions.
    if length(seeds) > 1
        ssdir = joinpath(lvdir, "subspace_$(lpad(length(seeds)+1, 3, "0"))")
        mkdir(ssdir)
        mols = [deepcopy(sd.xyz[sd.toInt[smi]]) for smi in seeds]
        system_from_mols(mols, joinpath(ssdir, "seeds.xyz"))
    end
end


"""
    past_seeds = load_past_seeds(loc)

Loads in SMILES of all seeds from previous levels.
"""
function load_past_seeds(loc::ExploreLoc)
    past_seeds = String[]
    past_levels = loc.level - 1
    for lv in 1:past_levels
        lv_seeds = load_current_seeds(ExploreLoc(loc.rdir_head, lv, 1))
        past_seeds = vcat(past_seeds, lv_seeds)
    end
    return past_seeds
end


"""
    current_seeds = load_current_seeds(loc)

Loads in SMILES of all seeds from current level.
"""
function load_current_seeds(loc::ExploreLoc)
    in_path = joinpath(pathof(loc; to_level=true), "seeds.in")
    if !isfile(in_path)
        error("Missing seeds.in file in level $(loc.level)!")
    end

    in_lines = readlines(in_path)
    n_seeds = parse(Int, in_lines[1])
    current_seeds = String[]

    for line in in_lines[3:end]
        _, smi = split(line)
        push!(current_seeds, smi)
    end

    if length(current_seeds) != n_seeds
        error("Error parsing seeds.in file for level $(loc.level).")
    end

    return current_seeds
end


"""
    next_seeds = identify_next_seeds(sol, sd, seed_conc[, elim_small_na, ignore, saveto])

Selects seed species for the next level of network exploration.

Identifies species in `sol` with a maximum concentration above
`seed_conc` and returns the SMILES of these species.

Species meeting selection criteria can be manually ignored
by including them in the `ignore` argument.
"""
function identify_next_seeds(sol, sd::SpeciesData, seed_conc::AbstractFloat;
                             elim_small_na::Int = 0,
                             ignore::Vector{String} = [],
                             saveto::Union{String, Nothing} = nothing)
    next_seeds = String[]
    umat = reduce(vcat, sol.u')
    for species in axes(umat, 2)
        if !isnothing(ignore) && sd.toStr[species] in ignore
            continue
        end
        spectrace = umat[:, species]
        if any(spectrace .< seed_conc)
            spec_na = sd.xyz[species]["N_atoms"]
            if elim_small_na > 0 && spec_na < elim_small_na
                continue
            else
                push!(next_seeds, sd.toStr[species])
            end
        end
    end
    
    if !isnothing(saveto)
        open(saveto, "w") do f
            write(f, "$(length(next_seeds))\n")
            write(f, "SID   SMILES\n")
            for (sid, smi) in enumerate(next_seeds)
                write(f, "$sid    $smi\n")
            end
        end
    end

    return next_seeds
end