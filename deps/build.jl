using Conda

function flush_log(logger=Base.global_logger())
    flush(logger.stream)
end

PROJECT_DIR = dirname(dirname(@__FILE__))

# Check git is installed.
gitcheck = run(Cmd(`which git`, ignorestatus=true))
if gitcheck.exitcode == 1
    error("Git not found on the system. Install Git before running this script.")
end

# Get dependencies
depsdir = joinpath(PROJECT_DIR, "submodules")
cdedir = joinpath(depsdir, "cde")
obcrdir = joinpath(depsdir, "OBCanonicalRadicals")
if readdir(cdedir) != String[] && readdir(obcrdir) != String[]
    @info string("Submodules already initialised, checking for updates...")
    gitcmd = run(Cmd(`git submodule update`, dir=PROJECT_DIR))
    if gitcmd.exitcode == 1
        error("Error updating Git submodules.")
    end
    @info string("Submodules up to date.")
else
    @info string("Obtaining and initialising submodules...")
    gitcmd = run(Cmd(`git submodule init`, dir=PROJECT_DIR))
    if gitcmd.exitcode == 1
        error("Error initialising Git submodules.")
    end
    gitcmd = run(Cmd(`git submodule update --init`, dir=PROJECT_DIR))
    if gitcmd.exitcode == 1
        error("Error initialising Git submodules.")
    end
    @info string("All submodules initialised.")
end
flush_log()

# Compile CDE.
@info string("CDE may require some configuration before compilation.")
cdemake = lowercase(get(ENV, "KINETICA_BUILD_CDE", "true"))
if cdemake == "true"
    @info string("Compiling CDE...")
    flush_log()
    cdemakecmd = run(Cmd(`make`, dir=cdedir))
    if cdemakecmd.exitcode == 1
        error("CDE failed to compile!")
    end
    @info string("CDE compilation complete.")
elseif cdemake == "false"
    @info string("ENV[\"KINETICA_BUILD_CDE\"] is false, skipping CDE compilation.")
else
    error("Unknown value for ENV[\"KINETICA_BUILD_CDE\"]: $(cdemake)")
end
flush_log()

# Ensure the Conda environment has the necessary dependencies.
@info "Setting up Python dependencies..."
flush_log()
Conda.add(["numpy", "openbabel"]; channel="conda-forge")
Conda.pip_interop(true)
Conda.pip("install", ["--no-deps", "-e", "$(obcrdir)"])
@info "Python setup complete."

@info string("KineticaCore setup complete.")