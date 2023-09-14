"""
    logger = start_log(logdir[, label])

Creates a `SimpleLogger` for logging to a file in `logdir`.

Filename of the log can be customised using `label`. All
new logfiles will have dates attached.

# Example
```julia
# Creates a new logfile in `logdir/MyLog_yymmdd-HHMMSS.log`
logger = start_log(logdir; label="MyLog")
```
"""
function start_log(logdir::String; label::String="Kinetica")
    if !isdir(logdir) mkdir(logdir) end
    date = Dates.format(now(), "yymmdd-HHMMSS")
    logio = open(joinpath(logdir, "$(label)_$(date).log"), "w")
    logger = SimpleLogger(logio)
    return logger
end

"""
    end_log(logger)

Closes the `IOStream` attached to `logger`.
"""
function end_log(logger::AbstractLogger)
    close(logger.stream)
end


"""
    flush(logger)

Flushes the `IOStream` attached to the `logger`.
"""
function Base.flush(logger::AbstractLogger)
    Base.flush(logger.stream)
end

"""
    flush_log()

Flushes the `IOStream` attached to the currently scoped logger.
"""
function flush_log()
    Base.flush(current_logger().stream)
end

"""
    with_global_logger(function)

Execute `function`, directing all log messages to the global logger.

Useful for temporarily bypassing task-local loggers when the
global logger is required, e.g. within progress bars.

# Example

```julia
foo(x) = @info "x = \$x"

foo("local") # Prints to closest locally-scoped logger

with_global_logger() do
    foo("global") # Prints to global logger
end
```
"""
function with_global_logger(@nospecialize(f::Function))
    with_logger(f, global_logger())
end