"""
    make_rcount(path::String)

Determines the number of reaction mechanisms generated.

Reads the number of generated reactions from the file at
`path`. If this file does not exist, creates it and sets
this number to zero. Returns this number.
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
    make_inert_file(path::String, inert_species::Vector{String})

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
    import_mechanism(loc::ExploreLoc, rcount[, max_molecularity=2, duplicate_reverse=true, unique_rxns=true])

Create a CRN's initial `SpeciesData` and `RxData` from a CDE generated mechanism(s).

Reads in the results of a CDE run at the `rcount` reaction
directory specified by `loc`. Returns a new `SpeciesData` and `RxData`
containing the unique species and reactions within these results,
provided these reactions do not exceed the maximum molecularity
set by `max_molecularity`, which defaults to only accepting
unimolecular and bimolecular reactions.

`duplicate_reverse` specifies whether to add the reverse reactions
of each of the imported reactions. `unique_rxns` specifies whether
to only add unique reactions to `RxData`.
"""
function import_mechanism(loc::ExploreLoc, rcount; 
        max_molecularity=2, duplicate_reverse=true, unique_rxns=true)
    rdir = pathof(loc)
    rsmis, rxyzs, rsys, psmis, pxyzs, psys, dHs = ingest_cde_run(rdir, rcount; 
                                                                 duplicate_reverse=duplicate_reverse)
    all_smis = vcat(reduce(vcat, rsmis), reduce(vcat, psmis))
    all_xyzs = vcat(reduce(vcat, rxyzs), reduce(vcat, pxyzs))
    sd = SpeciesData(all_smis, all_xyzs, loc.level)
    rd = RxData(sd, rsmis, psmis, rsys, psys, dHs, loc.level; 
                max_molecularity=max_molecularity, unique_rxns=unique_rxns)
    return sd, rd
end

"""
    import_mechanism!(sd::SpeciesData, rd::RxData, loc::ExploreLoc, rcount[, max_molecularity=2, duplicate_reverse=true, unique_rxns=true])

Extend a CRN's `SpeciesData` and `RxData` from a CDE generated mechanism(s).

Reads in the results of a CDE run at the `rcount` reaction
directory specified by `loc`. Extends `sd` and `rd` with the unique 
species and reactions within these results, provided these 
reactions do not exceed the maximum molecularity set by 
`max_molecularity`, which defaults to only accepting unimolecular
and bimolecular reactions.

`duplicate_reverse` specifies whether to add the reverse reactions
of each of the imported reactions. `unique_rxns` specifies whether
to only add unique reactions to `RxData`.
"""
function import_mechanism!(sd::SpeciesData, rd::RxData, loc::ExploreLoc, rcount;
        max_molecularity=2, duplicate_reverse=true, unique_rxns=true)
    rdir = pathof(loc)
    rsmis, rxyzs, rsys, psmis, pxyzs, psys, dHs = ingest_cde_run(rdir, rcount; 
                                                                 duplicate_reverse=duplicate_reverse)
    all_smis = vcat(reduce(vcat, rsmis), reduce(vcat, psmis))
    all_xyzs = vcat(reduce(vcat, rxyzs), reduce(vcat, pxyzs))
    push_unique!(sd, all_smis, all_xyzs, loc.level)
    push!(rd, sd, rsmis, psmis, rsys, psys, dHs, loc.level;
          max_molecularity=max_molecularity, unique_rxns=unique_rxns)
    return
end


"""
    import_network(rdir_head::String)

Imports a network from a level tree within `rdir_head`.

Recurses through level directories in `rdir_head`, then recurses
through subspace directories in each level. Within each subspace,
imports all mechanisms within each CDE run.

Returns the resulting network (an instance of `SpeciesData` and 
`RxData`).
"""
function import_network(rdir_head::String)
    @info "Importing all reactions in level tree under $(rdir_head)"; flush_log()
    level_dirs = readdir(rdir_head)
    level_dirs = level_dirs[startswith.(level_dirs, "level_")]
    n_levels = length(level_dirs)
    if length(n_levels) == 0
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
        xyz = frame_from_smiles(spec)
        push_unique!(sd, spec, xyz, 0)
    end

    # Loop through each level, adding each subspace.
    loc = ExploreLoc(rdir_head, 1, 1)
    for i in 1:n_levels
        reset_subspace!(loc)
        ss_dirs = readdir(pathof(loc; to_level=true))
        ss_dirs = ss_dirs[startswith.(ss_dirs, "subspace_")]
        n_subspaces = length(ss_dirs)
        for j in 1:n_subspaces
            rcount = make_rcount(joinpath(pathof(loc), "rcount"))
            for reac in 1:rcount
                import_mechanism!(sd, rd, loc, reac)
            end
            inc_subspace!(loc)
        end
        inc_level!(loc)
    end

    @info "Finished network import."
    @info "Network contains $(sd.n) species over $(rd.nr) reactions, explored over $(n_levels) levels.\n"
    flush_log()

    return sd, rd
end


"""
    cleanup_network(rdir_head::String)

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
    return
end


"""
    setup_level(loc::ExploreLoc, sd::SpeciesData, seeds::Vector{String})

Sets up a new level of iterative network exploration.

Uses the seed species in `seeds` to create a seed mapping
file `seeds.in` for the current level, denoted by `loc`. 
Generates a `seeds.xyz` for each subspace to be explores
within the level, including same-species and cross-species
subspaces.
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
    return
end


"""
    load_past_seeds(loc::ExploreLoc)

Loads in SMILES of all seeds from previous levels.

Returns an array of these SMILES.
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
    load_current_seeds(loc::ExploreLoc)

Loads in SMILES of all seeds from current level.

Returns an array of these SMILES.
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
    identify_next_seeds(sol, sd::SpeciesData, seed_conc<:AbstractFloat[, elim_small_na::Int=0, 
                        ignore::Vector{String}=[], saveto::Union{String, Nothing}=nothing])
    identify_next_seeds(sol, sd::SpeciesData[, elim_small_na::Int=0, ignore::Vector{String}=[],
                        saveto::Union{String, Nothing}=nothing])

Selects seed species for the next level of network exploration.

Identifies species in `sol` with a maximum concentration above
`seed_conc` and returns an array of the SMILES of these species.
If `seed_conc` is not provided, assumes all species in `sd` should
become seeds.

Species meeting selection criteria can be manually ignored
by including their SMILES in the `ignore` argument. Similarly,
if species below a certain number of atoms are not desired as
seeds (e.g. if they are too small to break down), this number
of atoms can be set with `elim_small_na`.

If `saveto` is set to a file path, the selected seeds and their
maximum concentrations are output to this file (usually
generated automatically as a 'seeds.out' file in the current 
level of exploration).
"""
function identify_next_seeds(sol, sd::SpeciesData, seed_conc::AbstractFloat;
                             elim_small_na::Int = 0,
                             ignore::Vector{String} = [],
                             saveto::Union{String, Nothing} = nothing)
    next_seeds = String[]
    next_seed_concs = Float64[]
    umat = reduce(vcat, sol.u')
    for species in axes(umat, 2)
        if !isnothing(ignore) && sd.toStr[species] in ignore
            continue
        end
        spectrace = umat[:, species]
        spec_max_conc = maximum(spectrace)
        if spec_max_conc >= seed_conc
            spec_na = sd.xyz[species]["N_atoms"]
            if elim_small_na > 0 && spec_na < elim_small_na
                continue
            else
                push!(next_seeds, sd.toStr[species])
                push!(next_seed_concs, spec_max_conc)
            end
        end
    end
    
    if !isnothing(saveto)
        max_smi_length = maximum(length.(next_seeds))
        open(saveto, "w") do f
            write(f, "$(length(next_seeds))\n")
            write(f, "SID   $(rpad("SMILES", max_smi_length))   Max. Conc.\n")
            for (sid, (smi, conc)) in enumerate(zip(next_seeds, next_seed_concs))
                write(f, "$(rpad(sid, 5)) $(rpad(smi, max_smi_length))   $conc\n")
            end
        end
    end

    return next_seeds
end

function identify_next_seeds(sol, sd::SpeciesData;
                             elim_small_na::Int = 0,
                             ignore::Vector{String} = [],
                             saveto::Union{String, Nothing} = nothing)
    next_seeds = String[]
    next_seed_concs = Float64[]
    umat = reduce(vcat, sol.u')
    for species in axes(umat, 2)
        if !isnothing(ignore) && sd.toStr[species] in ignore
            continue
        end
        spectrace = umat[:, species]
        spec_max_conc = maximum(spectrace)
        spec_na = sd.xyz[species]["N_atoms"]
        if elim_small_na > 0 && spec_na < elim_small_na
            continue
        else
            push!(next_seeds, sd.toStr[species])
            push!(next_seed_concs, spec_max_conc)
        end
    end
    
    if !isnothing(saveto)
        max_smi_length = maximum(length.(next_seeds))
        open(saveto, "w") do f
            write(f, "$(length(next_seeds))\n")
            write(f, "SID   $(rpad("SMILES", max_smi_length))   Max. Conc.\n")
            for (sid, (smi, conc)) in enumerate(zip(next_seeds, next_seed_concs))
                write(f, "$(rpad(sid, 5)) $(rpad(smi, max_smi_length))   $conc\n")
            end
        end
    end

    return next_seeds
end