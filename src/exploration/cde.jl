"""
    CDE(template_dir::String [, env_threads, cde_exec, sampling_seed, radius, nrxn, parallel_runs, parallel_exes, write_stdout, write_stderr, allow_errors])

CDE runner. Initialised through keyword-based struct, run through functor.

Struct contains fields for:

* CDE template directory (`template_dir`)
* Environmental multithreading number of threads (Optional, defaults to 1 thread; `env_threads`)
* Path to CDE executable (Optional, defaults to CDE packaged within CDE_jll; `cde_exec`)
* Seed for CDE's RNG (Optional, setting to 0 indicates seed should be random; `sampling_seed`)
* Radius for exploration of breakdown space (Optional, Default = 50; `radius`)
* Number of mechanisms to generate within a single CDE run (Optional, Default = 1; `nrxn`)
* Number of parallel CDE runs to execute (Optional, Default = 1; `parallel_runs`)
* Maximum number of parallel CDE executables to run at any time (Optional, Default = 1; `parallel_exes`)
* Whether to write CDE's stdout to file (Optional, Default = `false`; `write_stdout`)
* Whether to write CDE's stderr to file (Optional, Default = `false`; `write_stderr`)
* Whether to allow functions to continue running if CDE errors are detected (Optional, Default = `false`; `allow_errors`)

Additionally, some fields are usually modified within Kinetica,
and are not intended to be changed by users.

* Main reaction directory (`rdir`)
* XYZ file of starting molecule/material (`init_xyz`)
"""
@kwdef mutable struct CDE
    template_dir::String
    
    # Optional fields
    env_threads::Int = 1
    cde_exec::Union{String, Cmd} = run_cde()
    sampling_seed::Int = 0
    radius::Int = 50
    nrxn::Int = 1
    parallel_runs::Int = 1
    parallel_exes::Int = parallel_runs
    write_stdout::Bool = true
    write_stderr::Bool = false
    allow_errors::Bool = false

    # Fields usually handled by Kinetica
    rdir::String = "CHANGEME"
    init_xyz::String = "seeds.xyz"
end


"""
    (self::CDE)(rcount::Int)

Runs CDE for reaction `rcount` in directory `self.rdir`.

Serial runner, only spawns 1 CDE process.
"""
function (self::CDE)(rcount::Int)

    @info "--- Reaction $rcount ---"
    @info " - Starting new reaction mechanism generation."
    flush_log()

    # Prepare new directory structure.
    rxdir = joinpath(self.rdir, "reac_$(lpad(rcount, 5, "0"))")
    cp(self.template_dir, rxdir)
    cp(self.init_xyz, joinpath(rxdir, "Start.xyz"))

    # Prepare sampling input.
    self.sampling_seed == 0 ? seed = rand(1:100000) : seed = self.sampling_seed + rcount
    open(joinpath(rxdir, "input"), "a") do f
        write(f, "nmcrxn $(self.nrxn)\n")
        write(f, "nrxn $(self.radius)\n")
        write(f, "ranseed $seed\n")
    end
    @info "   - Set up single-ended breakdown path sampling in $rxdir."
    @info "   - Running sampling..."
    flush_log()

    # Run single-ended sampling for breakdown products.
    outfile = self.write_stdout ? joinpath(rxdir, "cde.out") : devnull
    errfile = self.write_stderr ? joinpath(rxdir, "cde.err") : devnull
    if self.cde_exec isa Cmd
        env_multithread!(self.cde_exec, self.env_threads)
        run(pipeline(Cmd(`$(self.cde_exec) input`, dir=rxdir), stdout=outfile, stderr=errfile))
    else
        mt_env = env_multithread(self.env_threads)
        run(pipeline(Cmd(`$(self.cde_exec) input`, env=mt_env, dir=rxdir), stdout=outfile, stderr=errfile))
    end

    # Check calculation ran correctly.
    success = true
    f = open(joinpath(rxdir, "input.log"), "r")
    outlines = readlines(f)
    for line in outlines
        if occursin("ERROR", line)
            @warn "Error in CDE run, check logs for more information"
            success = false
        end
    end
    close(f)
    if !ispath(joinpath(rxdir, "rxn_0001_step_0001.xyz"))
        @warn "Error in CDE run, no reaction steps found ($(rxdir))"
        success = false
    end

    # Report success and update reaction counter.
    if success
        @info "   - Sampling completed successfully!"
        open(joinpath(self.rdir, "rcount"), "w") do f
            write(f, "$(lpad(rcount, 5, "0"))")
        end
        flush_log()
        return true
    else
        if !self.allow_errors
            emsg = "Forbidden error in CDE run, stopping exploration."
            @error emsg
            throw(ErrorException(emsg))
        end
        @info "   - Sampling failed, removing directory."
        rm(rxdir; recursive=true)
        flush_log()
        return false
    end
end


"""
    (self::CDE)(rcountrange<:AbstractUnitRange)

Runs CDE for reactions in `rcountrange` in directory `cde.rdir`.

Can be used to run CDE calculations in serial or parallel, depending
on the values of `cde.parallel_procs` and `cde.parallel_exes`. if
these are greater than 1, will automatically split `length(rcountrange)`
calculations and run in parallel. Otherwise, runs calculations in serial. 
"""
function (self::CDE)(rcountrange::AbstractUnitRange)

    @info "--- Reactions $(rcountrange.start) - $(rcountrange.stop) ---"
    @info " - Starting new reaction mechanism generation."
    flush_log()

    # Prepare new directory structure.
    rxdirs = String[]
    for rc in rcountrange
        rxdir = joinpath(self.rdir, "reac_$(lpad(rc, 5, "0"))")
        push!(rxdirs, rxdir)
        cp(self.template_dir, rxdir)
        cp(self.init_xyz, joinpath(rxdir, "Start.xyz"))

        # Prepare sampling input.
        self.sampling_seed == 0 ? seed = rand(1:100000) : seed = self.sampling_seed + rc
        open(joinpath(rxdir, "input"), "a") do f
            write(f, "nmcrxn $(self.nrxn)\n")
            write(f, "nrxn $(self.radius)\n")
            write(f, "ranseed $seed\n")
        end
    end
    @info "   - Set up single-ended breakdown path sampling in $(length(rcountrange)) directories."
    @info "   - Running sampling..."
    flush_log()

    # Run single-ended sampling for breakdown products.
    parallel_procs = []
    if self.cde_exec isa Cmd
        env_multithread!(self.cde_exec, self.env_threads)
    else
        mt_env = env_multithread(self.env_threads)
    end
    for i in 1:length(rcountrange)
        outfile = self.write_stdout ? joinpath(rxdirs[i], "cde.out") : devnull
        errfile = self.write_stderr ? joinpath(rxdirs[i], "cde.err") : devnull
        if self.cde_exec isa Cmd
            push!(parallel_procs, pipeline(Cmd(`$(self.cde_exec) input`, dir=rxdirs[i]), 
                stdout=outfile, stderr=errfile))
        else
            push!(parallel_procs, pipeline(Cmd(`$(self.cde_exec) input`, env=mt_env, dir=rxdirs[i]), 
              stdout=outfile, stderr=errfile))
        end
    end

    parallel_run(parallel_procs; ntasks=self.parallel_exes)

     # Check calculations ran correctly.
    success = [true for _ in 1:length(rcountrange)]
    for i in 1:length(rcountrange)
        f = open(joinpath(rxdirs[i], "input.log"), "r")
        outlines = readlines(f)
        for line in outlines
            if occursin("ERROR", line) 
                success[i] = false 
                @warn "Error in CDE run $i, check logs for more information."
            end
        end
        close(f)
        if !ispath(joinpath(rxdirs[i], "rxn_0001_step_0001.xyz"))
            success[i] = false
            @warn "Error in CDE run $i, no reaction steps found ($(rxdirs[i]))"
        end
    end

    # Report success and update reaction counter.
    if all(success)
        @info "   - Sampling completed successfully!"
        open(joinpath(self.rdir, "rcount"), "w") do f
            write(f, "$(lpad(rcountrange.stop, 5, "0"))")
        end
    else
        if !self.allow_errors
            emsg = "Forbidden error in at least one CDE run, stopping exploration."
            @error emsg
            throw(ErrorException(emsg))
        end
        for (i, s) in enumerate(success)
            if !s
                @info " - Sampling failed in CDE run $i, removing directory."
                rm(rxdirs[i]; recursive=true)
            end
        end
        @info "   - Reordering reactions."
        counter = 0
        for i in 1:length(rcountrange)
            if success[i]
                counter += 1
                i != counter && mv(rxdirs[i], rxdirs[counter])
            end
        end
        open(joinpath(self.rdir, "rcount"), "w") do f
            write(f, "$(lpad(rcountrange.start+count(success)-1, 5, "0"))")
        end
    end
        
    flush_log()
    return rcountrange.start+count(success)-1

end


"""
    ingest_cde_run(rdir::String, rcount[, fix_radicals=true])

Reads in the results from a CDE run.

Separates out fragment species from each available reaction's
reactants and products, forming arrays of their SMILES strings
and ExtXYZ geometries.

OBCanonicalRadicals can be enabled to tidy up radical SMILES
using the `fix_radicals` parameter.

Returns `reac_smis, reac_xyzs, reac_systems, prod_smis, prod_xyzs, prod_systems, dH`,
where:

* `reac_smis` and `prod_smis` are arrays of the SMILES of
each reaction's reactants and products;
* `reac_xyzs` and `prod_xyzs` are their corresponding 
geometries as ExtXYZ frames; 
* `reac_systems` and `prod_systems` are the ExtXYZ frames
of the systems of molecules that came out of CDE;
* `dH` is an array of reaction energies.
"""
function ingest_cde_run(rdir::String, rcount; fix_radicals=true)
    rxdir = joinpath(rdir, "reac_$(lpad(rcount, 5, "0"))")
    @debug "Reading in mechanism step xyz files."

    # Read in all reactions as 2-frame trajectories.
    rxfiles = readdir(rxdir)
    rxfiles = rxfiles[startswith.(rxfiles, "rxn_")]
    n_reacs = length(rxfiles)
    reacs = []
    prods = []
    dH = zeros(Float64, n_reacs)
    for i in 1:n_reacs
        reaction = read_frames(joinpath(rxdir, rxfiles[i]), 1:2)

        # Separate out reacs and prods, calculate dH.
        push!(reacs, reaction[1])
        push!(prods, reaction[2])
        dH[i] = reaction[2]["info"]["energy"] - reaction[1]["info"]["energy"]
    end

    @debug "Extracting fragment species from reactions."
    reac_smis = Vector{String}[]
    reac_xyzs = Vector{Dict{String, Any}}[]
    reac_systems = Dict{String, Any}[]
    for reac in reacs
        smis, xyzs = ingest_xyz_system(frame_to_xyz(reac); fix_radicals)
        push!(reac_smis, smis)
        push!(reac_xyzs, xyzs)
        push!(reac_systems, reac)
    end
    prod_smis = Vector{String}[]
    prod_xyzs = Vector{Dict{String, Any}}[]
    prod_systems = Dict{String, Any}[]
    for prod in prods
        smis, xyzs = ingest_xyz_system(frame_to_xyz(prod); fix_radicals)
        push!(prod_smis, smis)
        push!(prod_xyzs, xyzs)
        push!(prod_systems, prod)
    end

    # Add in all reverse reactions.
    reac_smis = vcat(reac_smis, prod_smis)
    prod_smis = vcat(prod_smis, reac_smis)
    dH = vcat(dH, -dH)
    reac_systems = vcat(reac_systems, prod_systems)
    prod_systems = vcat(prod_systems, reac_systems)
    @debug "Read in $(n_reacs*2) reactions."

    return reac_smis, reac_xyzs, reac_systems, prod_smis, prod_xyzs, prod_systems, dH
end