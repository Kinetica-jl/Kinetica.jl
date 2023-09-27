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
    @info "Importing all reactions in level tree under $(rdir_head)"
    level_dirs = readdir(rdir_head)
    level_dirs = level_dirs[startswith.(level_dirs, "level_")]
    if length(level_dirs) == 0
        error("ERROR: No network levels found in rdir_head.")
    end
    
    # Initialise network with blank sd and rd.
    sd = SpeciesData{Int64}(
        Dict{String, Int64}(), Dict{Int64, String}(),
        0, 
        Dict{Int64, Dict{String, Any}}(), Dict()
    )
    rd = RxData{Int64, Float64}(
        0, 
        Vector{String}[], Vector{String}[],
        Vector{Int64}[], Vector{Int64}[], 
        Vector{Int64}[], Vector{Int64}[], 
        Float64[], Vector{UInt8}[]
    )

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
    @info "Network contains $(sd.n) species over $(rd.nr) reactions."
    flush_log()

    return sd, rd
end