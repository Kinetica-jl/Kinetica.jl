"""
    cde = CDE(rdir, tdir, init_xyz, env[, sampling_seed, radius])

CDE runner. Initialised through struct, run through functor.

Struct contains fields for:
* Main reaction directory (`rdir`)
* CDE template directory (`tdir`)
* XYZ file of starting molecule/material (`init_xyz`)
* Environment variable dictionary, modified for running CDE calculations (`env`)
* Path to CDE executable (Optional, defaults to CDE compiled by KineticaCore's build process; `cde_exec`)
* Seed for CDE's RNG (Optional; `sampling_seed`)
* Radius for exploration of breakdown space (Optional, Default = 50; `radius`)
* Number of mechanisms to generate within a single CDE run (Optional, Default = 1; `nrxn`)
* Number of parallel CDE runs to execute (Optional, Default = 1; `parallel_runs`)
* Maximum number of parallel CDE executables to run at any time (Optional, Default = 1; `parallel_exes`)
* Whether to write CDE's stdout to file (Optional, Default = `false`; `write_stdout`)
* Whether to write CDE's stderr to file (Optional, Default = `false`; `write_stderr`)
* Whether to allow functions to continue running if CDE errors are detected (Optional, Default = `false`; `allow_errors`)

--- Functor definition ---

    cderun(rcount)

Functor for running CDE with the parameters set in the struct.

`rcount` should the the next available reaction ID, not the current
one in the rcount file.

Sets up directories for CDE calculation of the next available mechanism,
then runs the calculation and reports any errors.

"""
Base.@kwdef mutable struct CDE
    # Mandatory fields.
    rdir::String
    tdir::String
    init_xyz::String
    env::Dict

    # Optional fields.
    cde_exec::String = joinpath(dirname(dirname(@__FILE__)), "submodules/cde/bin/cde.x")
    sampling_seed::Int = 0
    radius::Int = 50
    nrxn::Int = 1
    parallel_runs::Int = 1
    parallel_exes::Int = 1
    write_stdout::Bool = true
    write_stderr::Bool = false
    allow_errors::Bool = false
end


"""
    cde(rcount)

Runs CDE for reaction `rcount` in directory `cde.rdir`.

Serial runner, only spawns 1 CDE process.
"""
function (self::CDE)(rcount::Int)

    @info "--- Reaction $rcount ---"
    @info "Starting new reaction mechanism generation."
    flush_log()

    # Prepare new directory structure.
    rxdir = joinpath(self.rdir, "reac_$(lpad(rcount, 5, "0"))")
    cp(self.tdir, rxdir)
    cp(self.init_xyz, joinpath(rxdir, "Start.xyz"))

    # Prepare sampling input.
    self.sampling_seed == 0 ? seed = rand(1:100000) : seed = self.sampling_seed
    open(joinpath(rxdir, "input"), "a") do f
        write(f, "nmcrxn $(self.nrxn)\n")
        write(f, "nrxn $(self.radius)\n")
        write(f, "ranseed $seed\n")
    end
    @info "Set up single-ended breakdown path sampling in $rxdir."
    @info "Running sampling..."
    flush_log()

    # Run single-ended sampling for breakdown products.
    outfile = self.write_stdout ? joinpath(rxdir, "cde.out") : devnull
    errfile = self.write_stderr ? joinpath(rxdir, "cde.err") : devnull
    run(pipeline(Cmd(`$(self.cde_exec) input`, env=self.env, dir=rxdir), stdout=outfile, stderr=errfile))

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
        @info "Sampling completed successfully!\n"
        open(joinpath(self.rdir, "rcount"), "w") do f
            write(f, "$(lpad(rcount, 5, "0"))")
        end
        flush_log()
        return true
    else
        if !self.allow_errors
            emsg = "Forbidden error in CDE run, stopping exploration."
            @error emsg
            throw(ErrorException("Forbidden error in CDE run, stopping exploration."))
        end
        @info "Sampling failed, removing directory.\n"
        rm(rxdir; recursive=true)
        flush_log()
        return false
    end
end


"""
    cde(rcountrange)

Runs CDE for reactions in `rcountrange` in directory `cde.rdir`.

Can be used to run CDE calculations in serial or parallel, depending
on the values of `cde.parallel_procs` and `cde.parallel_exes`. if
these are greater than 1, will automatically split `length(rcountrange)`
calculations and run in parallel. Otherwise, runs calculations in serial. 
"""
function (self::CDE)(rcountrange::AbstractUnitRange)

    @info "--- Reactions $(rcountrange.start) - $(rcountrange.stop) ---"
    @info "Starting new reaction mechanism generation."
    flush_log()

    # Prepare new directory structure.
    rxdirs = String[]
    for rc in rcountrange
        rxdir = joinpath(self.rdir, "reac_$(lpad(rc, 5, "0"))")
        push!(rxdirs, rxdir)
        cp(self.tdir, rxdir)
        cp(self.init_xyz, joinpath(rxdir, "Start.xyz"))

        # Prepare sampling input.
        self.sampling_seed == 0 ? seed = rand(1:100000) : seed = self.sampling_seed
        open(joinpath(rxdir, "input"), "a") do f
            write(f, "nmcrxn $(self.nrxn)\n")
            write(f, "nrxn $(self.radius)\n")
            write(f, "ranseed $seed\n")
        end
    end
    @info "Set up single-ended breakdown path sampling in $(length(rcountrange)) directories."
    @info "Running sampling..."
    flush_log()

    # Run single-ended sampling for breakdown products.
    parallel_procs = []
    for i in 1:length(rcountrange)
        outfile = self.write_stdout ? joinpath(rxdirs[i], "cde.out") : devnull
        errfile = self.write_stderr ? joinpath(rxdirs[i], "cde.err") : devnull
        push!(parallel_procs, pipeline(Cmd(`$(self.cde_exec) input`, env=self.env, dir=rxdirs[i]), 
              stdout=outfile, stderr=errfile))
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
        @info "Sampling completed successfully!\n"
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
                @info "Sampling failed in CDE run $i, removing directory."
                rm(rxdirs[i]; recursive=true)
            end
        end
        @info "Reordering reactions.\n"
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

