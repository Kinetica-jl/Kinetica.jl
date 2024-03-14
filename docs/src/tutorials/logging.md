# Logging

Kinetica makes use of multiple logging functions from the [JuliaLogging](https://julialogging.github.io) organisation to handle writing logs of varying detail levels to the console and to file, as well as to handle various optional progress bars for tracking ODE solution progress.

By default, running Kinetica functions such as [`explore_network`](@ref) will output a text stream containing important information about the task at hand to `stdout`. An example of this can be seen in the ['Running the Simulation' section of Getting Started](@ref "Running the Simulation"). This will always log at the `Info` level, so detailed debugging messages will not be present.

To log to a file and/or enable `Debug`-level logging, Kinetica provides a shorthand function for setting up the correct logger, [`start_log`](@ref):

```julia
using Logging: with_logger, Debug, Info
logger = start_log("./"; min_level=Info, label="MyLog")
```

This creates a nicely formatted logger of the requested level at the path `./MyLog_yymmdd-HHMMSS.log`, inserting the date and time of creation into the `yymmdd` and `HHMMSS` fields respectively. If debug logging is required, `min_level=Debug` should be set instead.

To use this logger on a set of expressions, correctly formatting log messages and sending them to the requested log file, it needs to wrap the expressions in a [`with_logger`](https://julialogging.github.io/reference/logging/#Logging.with_logger) function:

```julia
with_logger(logger) do
    global res = explore_network(exploremethod, solvemethod, "./my_CRN_out")
end
```

Note the use of the `global` keyword here - if this is a top-level script, the [`with_logger`](https://julialogging.github.io/reference/logging/#Logging.with_logger) function will create a new scope that contains the variable `res`, which will be lost to the global (script-level) scope once it is exited. By making this variable `global`, it ensures we don't lose simulation results while logging!

Once you've finished logging within a script, call [`end_log`](@ref) to safely close the log file:

```julia
end_log(logger)
```

## Progress Bars

When solving ODEs, Kinetica can both make use of [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)s native progress bar implementation and set up its own, depending on the solve type requested. Both depend on [TerminalLoggers.jl](https://julialogging.github.io/reference/terminalloggers) in the background, but if we were to wrap a progress bar into a text log, we'd get a bit of a mess!

Kinetica works around this by always passing progress bars within ODE solution calls to the *global logger*. While the `logger` variable we defined above is a *local logger* that needs to be called within a [`with_logger`](https://julialogging.github.io/reference/logging/#Logging.with_logger) wrapper, the global logger is specified over the entire session and can be accessed from anywhere in the code. It is therefore enough to define the global logger at the start of a Kinetica script, and the solvers will handle the rest. As in [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), this is done by specifying the following:

```julia
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
```

To enable progress bars during ODE solution, the `progress` parameter of [`ODESimulationParams`](@ref) needs to be set to `true`.