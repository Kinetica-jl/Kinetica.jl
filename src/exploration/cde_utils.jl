"""
    parallel_run(commands[, ntasks])

Asynchronously runs `length(commands)` shell commands over `ntasks` processes.
"""
function parallel_run(commands; ntasks::Int=1)
    request = Channel() do request
        for cmd in commands
            put!(request, cmd)
        end
    end
    @sync for _ in 1:ntasks
        @async try
            foreach(run, request)
        finally
            close(request)  # shutdown on error
        end
    end
end


"""
    env = env_multithread(nthreads)

Set up environment variables for CDE calculations.

Sets OMP/MKL thread counts (for xTB calculations within CDE).
"""
function env_multithread(nthreads::Int)
    env = copy(ENV)

    env["OMP_NUM_THREADS"] = "$nthreads"
    env["MKL_NUM_THREADS"] = "$nthreads"
    env["MKL_DYNAMIC"] = "FALSE"

    return env
end

"""
    env_multithread!(cmd, nthreads)

Set up environment variables for CDE calculations by modifying an existing Cmd.

Sets OMP/MKL thread counts (for xTB calculations within CDE).
"""
function env_multithread!(cmd::Cmd, nthreads::Int)
    push!(cmd.env, "OMP_NUM_THREADS=$nthreads")
    push!(cmd.env, "MKL_NUM_THREADS=$nthreads")
    push!(cmd.env, "MKL_DYNAMIC=FALSE")
end