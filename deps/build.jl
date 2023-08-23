PROJECT_DIR = dirname(dirname(@__FILE__))

# Check git is installed.
gitcheck = run(Cmd(`which git`, ignorestatus=true))
if gitcheck.exitcode == 1
    error("Git not found on the system. Install Git before running this script.")
end

# Get dependencies
depsdir = joinpath(PROJECT_DIR, "submodules")
cdedir = joinpath(depsdir, "cde")
if readdir(cdedir) != String[]
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

# Compile CDE.
@info string("CDE may require some configuration before compilation.")
@info string("Compiling CDE...")
cdemakecmd = run(Cmd(`make`, dir=cdedir))
if cdemakecmd.exitcode == 1
    error("CDE failed to compile!")
end

@info string("CDE compilation complete.")
@info string("KineticaCore setup complete.")