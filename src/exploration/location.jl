"""
    loc = ExploreLoc(rdir_head, level, subspace)

Container for iterative exploration location data.
"""
mutable struct ExploreLoc
    rdir_head::String
    level::Int64
    subspace::Int64
end

inc_level!(loc::ExploreLoc) = loc.level += 1
inc_subspace!(loc::ExploreLoc) = loc.subspace += 1
dec_level!(loc::ExploreLoc) = loc.level -= 1
dec_subspace!(loc::ExploreLoc) = loc.subspace -= 1
reset_subspace!(loc::ExploreLoc) = loc.subspace = 1

"""
    path = pathof(loc[, to_level])

Returns the path of the current subspace in the current level in the head directory specified by `loc`.

If `to_level=true`, only returns the path to the current level,
excluding the current subspace.
"""
function Base.pathof(loc::ExploreLoc; to_level::Bool=false)
    if to_level
        return joinpath(loc.rdir_head, "level_$(lpad(loc.level, 3, "0"))")
    else
        return joinpath(loc.rdir_head, "level_$(lpad(loc.level, 3, "0"))", "subspace_$(lpad(loc.subspace, 3, "0"))")
    end
end


"""
"""
function find_current_loc(rdir_head::String)
    level_dirs = readdir(rdir_head)
    level_dirs = level_dirs[startswith.(level_dirs, "level_")]
    if length(level_dirs) == 0
        @info "No network levels found in $(rdir_head), starting network exploration from scratch."
        loc = ExploreLoc(rdir_head, 0, 1)
        return loc
    end

    curr_level_dir = level_dirs[end]
    level = parse(Int, split(curr_level_dir, "_")[end])
    if !isfile(joinpath(rdir_head, curr_level_dir, "seeds.in"))
        @info "No seeds.in found in level $level, continuing from previous level."
        curr_level_dir = level_dirs[end-1]
        level -= 1
    end
    rdir_level_dir = joinpath(rdir_head, curr_level_dir)

    subspace_dirs = readdir(rdir_level_dir)
    subspace_dirs = subspace_dirs[startswith.(subspace_dirs, "subspace_")]
    if length(subspace_dirs) == 0
        @info "No subspaces found in level $level, starting level exploration from scratch."
        loc = ExploreLoc(rdir_head, level, 1)
        return loc
    end

    subspace = 1
    for (i, subspace_dir) in enumerate(subspace_dirs)
        subspace = i
        curr_subspace_dir = joinpath(rdir_level_dir, subspace_dir)
        if !isfile(joinpath(curr_subspace_dir, "isconv"))
            @info "Current exploration location: Level $level, Subspace $subspace"
            loc = ExploreLoc(rdir_head, level, subspace)
            return loc
        end
    end

    @warn "All subspaces in level $level are converged!"
    @info "Current exploration location: Level $level, Subspace $subspace"
    loc = ExploreLoc(rdir_head, level, subspace)
    return loc
end