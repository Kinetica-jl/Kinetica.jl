abstract type AbstractSolveMethod end

abstract type AbstractODESolveMethod <: AbstractSolveMethod end
abstract type AbstractSSASolveMethod <: AbstractSolveMethod end

"""
Static kinetic CRN solver type.

Combines all parameter inputs, condtions and
a calculator into a single type to be passed
to the solver.

All conditions in the provided `ConditionSet`
must be static (defined as a single number),
and also must be compatible with the calculator.
"""
struct StaticODESolve <: AbstractODESolveMethod
    pars::ODESimulationParams
    conditions::ConditionSet
    calculator::AbstractKineticCalculator
    function StaticODESolve(pars, conditions, calculator)
        if !isstatic(conditions)
            throw(ArgumentError("All conditions must be static to run a StaticODESolve."))
        elseif !has_conditions(calculator, conditions.symbols)
            throw(ArgumentError("Calculator does not support all of the provided conditions."))
        else
            return new(pars, conditions, calculator)
        end
    end
end

"""
Variable kinetic CRN solver type.

Combines all parameter inputs, condtions and
a calculator into a single type to be passed
to the solver.

Conditions in the provided `ConditionSet` must
be compatible with the calculator and can be 
a combination of static and variable. However,
this will throw an error if all conditions are
static, as a `StaticODESolve` should be used
instead.
"""
struct VariableODESolve <: AbstractODESolveMethod
    pars::ODESimulationParams
    conditions::ConditionSet
    calculator::AbstractKineticCalculator
    function VariableODESolve(pars, conditions, calculator)
        if !has_conditions(calculator, conditions.symbols)
            throw(ArgumentError("Calculator does not support all of the provided conditions."))
        else
            return new(pars, conditions, calculator)
        end
    end
end




"""
    sol = solve_network(method::StaticODESolve, rd, species)

Solve a network with static kinetics.

Automatically dispatches to the correct method based
on the value of `method.pars.solve_chunks`, as chunkwise
solution requires a significantly different approach.
"""
function solve_network(method::StaticODESolve, rd::RxData, species::SpeciesData)
    split_method = method.pars.solve_chunks ? :chunkwise : :complete
    sol = solve_network(method, rd, species, Val(split_method))
    return sol
end

function solve_network(method::StaticODESolve, rd::RxData, species::SpeciesData, ::Val{:complete})
    if hasfield(method.calculator, :setup_network!)
        @info " - Setting up network for rate calculation"
        method.calculator.setup_network!(rd, species, method.pars)
    end

    @info " - Calculating rate constants"
    condition_map = [sym => profile.value for (sym, profile) in zip(method.calculator.symbols, method.calculator.profiles)]
    rates = method.calculator(; condition_map...)
    apply_low_k_cutoff!(rd, method.pars, rates)

    @info " - Setting up ReactionSystem"
    @parameters k[1:rd.nr]
    @variables t (spec(t))[1:species.n]

    u0 = make_u0(species, method.pars)
    u0map = Pair.(collect(spec), u0)
    pmap = Pair.(collect(k), rates)

    rs = make_rs(k, spec, t, rd)

    @info " - Formulating ODEProblem"
    @info "   - Sparse? $(method.pars.sparse)"
    @info "   - Analytic Jacobian? $(method.pars.jac)"
    oprob = ODEProblem(rs, u0map, method.pars.tspan, pmap;
        jac=method.pars.jac, sparse=method.pars.sparse)
    solvecall_kwargs = Dict{Symbol, Any}(
        :progress => method.opars.progress,
        :progress_steps => 10,
        :abstol => method.pars.abstol,
        :reltol => method.pars.reltol,
        :dtmin => eps(eltype(method.pars.tspan)),
        :maxiters => method.pars.maxiters,
        :saveat => isnothing(method.pars.save_interval) ? eltype(method.pars.tspan)[] : method.pars.save_interval,
        :kwargshandle => KeywordArgError
    )
    if method.pars.ban_negatives
        solvecall_kwargs[:isoutofdomain] = (u,p,t)->any(x->x<0,u)
    end

    integ = init(oprob, method.pars.solver(; method.pars.solver_kwargs...); solvecall_kwargs...)
    sol = adaptive_solve!(integ, method.pars, solvecall_kwargs; print_status=true)

    return sol
end

function solve_network(method::StaticODESolve, rd::RxData, species::SpeciesData, ::Val{:chunkwise})
    if hasfield(method.calculator, :setup_network!)
        @info " - Setting up network for rate calculation."
        method.calculator.setup_network!(rd, species)
    end

    @info " - Calculating rate constants."
    condition_map = [sym => profile.value for (sym, profile) in zip(method.calculator.symbols, method.calculator.profiles)]
    rates = method.calculator(; condition_map...)
    apply_low_k_cutoff!(rd, method.pars, rates)

    @info " - Setting up ReactionSystem"
    @parameters k[1:rd.nr]
    @variables t (spec(t))[1:species.n]

    u0 = make_u0(species, method.pars)
    u0map = Pair.(collect(spec), u0)
    pmap = Pair.(collect(k), rates)

    rs = make_rs(k, spec, t, rd)

    @info " - Formulating ODEProblem"
    @info "   - Sparse? $(method.pars.sparse)"
    @info "   - Analytic Jacobian? $(method.pars.jac)"
    tType = eltype(method.pars.tspan)
    uType = eltype(u0)
    local_tspan = (tType(0.0), tType(method.pars.solve_chunkstep))
    oprob = ODEProblem(rs, u0map, local_tspan, pmap;
        jac=method.pars.jac, sparse=method.pars.sparse)

    # Determine how many solution chunks will be required.
    n_chunks_reqd = Int(method.pars.tspan[2] / method.pars.solve_chunkstep)
    save_interval = isnothing(method.pars.save_interval) ? method.pars.solve_chunkstep : method.pars.save_interval
    saveat_local = collect(0.0:pars.save_interval:opars.solve_chunkstep)

    # Allocate final solution arrays.
    size_final = (length(saveat_local)-1)*n_chunks_reqd + 1
    u_final = zeros(uType, size_final, length(u0))
    t_final = zeros(tType, size_final)

    # Initialise integrator.
    solvecall_kwargs = Dict{Symbol, Any}(
        :progress => false,
        :abstol => method.pars.abstol,
        :reltol => method.pars.reltol,
        :dtmin => eps(method.pars.solve_chunkstep),
        :maxiters => method.pars.maxiters,
        :saveat => isnothing(method.pars.save_interval) ? eltype(method.pars.tspan)[] : method.pars.save_interval,
        :kwargshandle => KeywordArgError
    )
    if method.pars.ban_negatives
        solvecall_kwargs[:isoutofdomain] = (u,p,t)->any(x->x<0,u)
    end
    integ = init(oprob, method.pars.solver(; method.pars.solver_kwargs...); solvecall_kwargs...)

    # Set up progress bar (if required).
    if method.pars.progress
        pbar_sid = uuid4()
        with_global_logger() do
            @info Progress(pbar_sid, name="Chunkwise ODE")
        end
    end

    # Loop over the solution chunks needed to generate the full solution.
    for nc in 0:n_chunks_reqd-1
        # Reinitialise the integrator at the current concentrations.
        reinit!(integ, integ.sol.u[end])

        adaptive_solve!(integ, method.pars, solvecall_kwargs)
        if method.pars.progress 
            with_global_logger() do
                @info Progress(pbar_sid, (nc+1)/n_chunks_reqd)
            end
        end
        
        # Determine where to place the loop's results in the final solution arrays.
        start_idx = nc*(length(saveat_local)-1) + 1
        end_idx = start_idx + length(saveat_local) - 2

        # If on the last loop, include final timestep.
        umat = reduce(vcat, integ.sol.u')
        if nc == n_chunks_reqd-1
            end_idx += 1
            u_final[start_idx:end_idx, :] = umat[:, :]
            t_final[start_idx:end_idx] = integ.sol.t .+ (nc*method.pars.solve_chunkstep)
        # Otherwise, insert all but the final saved timesteps
        else
            u_final[start_idx:end_idx, :] = umat[1:end-1, :]
            t_final[start_idx:end_idx] = integ.sol.t[1:end-1] .+ (nc*method.pars.solve_chunkstep)
        end
    end

    if method.pars.progress 
        with_global_logger() do 
            @info Progress(pbar_sid, done=true)
        end 
    end

    # Somehow remake a normal solution object.

end