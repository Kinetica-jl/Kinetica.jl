"""
    cde = CDE(rdir, template_dir, init_xyz, env[, sampling_seed, radius])

CDE runner. Initialised through struct, run through functor.

Struct contains fields for:

* CDE template directory (`template_dir`)
* Environment variable dictionary, modified for running CDE calculations (Optional, defaults to current environment; `env`)
* Path to CDE executable (Optional, defaults to CDE compiled by KineticaCore's build process; `cde_exec`)
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
    env::Dict = ENV
    cde_exec::String = joinpath(dirname(dirname(dirname(@__FILE__))), "submodules/cde/bin/cde.x")
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
    cde(rcount)

Runs CDE for reaction `rcount` in directory `cde.rdir`.

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
    self.sampling_seed == 0 ? seed = rand(1:100000) : seed = self.sampling_seed
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
    cde(rcountrange)

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
        self.sampling_seed == 0 ? seed = rand(1:100000) : seed = self.sampling_seed
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
    reac_smis, reac_xyzs, prod_smis, prod_xyzs, dH = ingest_cde_run(rdir, rcount[, fix_radicals])

Reads in the results from a CDE run.

Separates out fragment species from each available reaction's
reactants and products, forming arrays of their SMILES strings
and ExtXYZ geometries.

OBCanonicalRadicals can be enabled to tidy up radical SMILES
using the `fix_radicals` parameter.
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
    tmp_xyz = joinpath(rxdir, "kinetica_tmp.xyz")

    reac_smis = Vector{String}[]
    reac_xyzs = []
    for reac in reacs
        write_frame(tmp_xyz, reac)
        smis, xyzs = ingest_xyz_system(tmp_xyz; fix_radicals)
        push!(reac_smis, smis)
        push!(reac_xyzs, xyzs)
    end

    prod_smis = Vector{String}[]
    prod_xyzs = []
    for prod in prods
        write_frame(tmp_xyz, prod)
        smis, xyzs = ingest_xyz_system(tmp_xyz; fix_radicals)
        push!(prod_smis, smis)
        push!(prod_xyzs, xyzs)
    end

    rm(tmp_xyz)

    # Add in all reverse reactions.
    reac_smis = vcat(reac_smis, prod_smis)
    prod_smis = vcat(prod_smis, reac_smis)
    dH = vcat(dH, -dH)
    @debug "Read in $(n_reacs*2) reactions."

    return reac_smis, reac_xyzs, prod_smis, prod_xyzs, dH
end