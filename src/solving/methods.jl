abstract type AbstractSolveMethod end

abstract type AbstractODESolveMethod <: AbstractSolveMethod end
abstract type AbstractSSASolveMethod <: AbstractSolveMethod end


struct StaticODESolve <: AbstractODESolveMethod
    pars::ODESimulationParams
    conditions::ConditionSet
    calculator::AbstractKineticCalculator
    filter::RxFilter
    function StaticODESolve(pars, conditions, calculator, filter)
        if !isstatic(conditions)
            throw(ArgumentError("All conditions must be static to run a StaticODESolve."))
        elseif !has_conditions(calculator, conditions.symbols)
            throw(ArgumentError("Calculator does not support all of the provided conditions."))
        else
            return new(pars, conditions, calculator, filter)
        end
    end
end

"""
    StaticODESolve(pars::ODESimulationParams, conditions::ConditionSet, calculator<:AbstractKineticCalculator[, filter::RxFilter])

Static kinetic CRN solver type.

Combines all parameter inputs, conditions, a calculator and
a set of reaction filters into a single type to be passed to
the solver.

All conditions in the provided `ConditionSet` must be static 
(defined as a single number), and also must be compatible with
the provided calculator.

Allows for optional specification of `RxFilter`. If not defined,
creates a filter that allows all reactions.
"""
function StaticODESolve(pars::ODESimulationParams, conditions::ConditionSet, calculator::AbstractKineticCalculator)
    return StaticODESolve(pars, conditions, calculator, RxFilter())
end


struct VariableODESolve <: AbstractODESolveMethod
    pars::ODESimulationParams
    conditions::ConditionSet
    calculator::AbstractKineticCalculator
    filter::RxFilter
    function VariableODESolve(pars, conditions, calculator, filter)
        if !has_conditions(calculator, conditions.symbols)
            throw(ArgumentError("Calculator does not support all of the provided conditions."))
        elseif !conditions.discrete_updates && !allows_continuous(calculator)
            throw(ArgumentError("Calculator does not support continuous rate updates in simulations."))
        else
            return new(pars, conditions, calculator, filter)
        end
    end
end

"""
    VariableODESolve(pars::ODESimulationParams, conditions::ConditionSet, calculator<:AbstractKineticCalculator[, filter::RxFilter])

Variable kinetic CRN solver type.

Combines all parameter inputs, conditions, a calculator and
a set of reaction filters into a single type to be passed to
the solver.

Conditions in the provided `ConditionSet` must be compatible
with the calculator and can be a combination of static and 
variable. However, this will throw an error if all conditions 
are static, as a `StaticODESolve` should be used instead.

Allows for optional specification of `RxFilter`. If not defined,
creates a filter that allows all reactions.
"""
function VariableODESolve(pars::ODESimulationParams, conditions::ConditionSet, calculator::AbstractKineticCalculator)
    return VariableODESolve(pars, conditions, calculator, RxFilter())
end


"""
    solve_network(method::StaticODESolve, sd::SpeciesData, rd::RxData[, copy_network=true, return_integrator=false])

Solve a network with static kinetics.

Automatically dispatches to the correct method based
on the value of `method.pars.solve_chunks`, as chunkwise
solution requires a significantly different approach.

Setting `copy_network=true` generates a `deepcopy` of
the original network in `rd` and `sd` and uses these in
the solution, to avoid side effects from calculators
modifying the original network that is passed in. The
copied (modified) network is returned as part of the 
resulting `ODESolveOutput`.

Setting `return_integrator=true` sets up and returns the
underlying integrator without solving (i.e. at 
`t = method.pars.tspan[1]`), allowing for manual stepping
through the solution. Note that chunkwise solutions implement
many reinitialisations of this integrator, which will have
to be mirrored in the calling script to get the same results.
"""
function solve_network(method::StaticODESolve, sd::SpeciesData, rd::RxData; 
                       copy_network=true, return_integrator=false)
    if copy_network
        sd_active = deepcopy(sd)
        rd_active = deepcopy(rd)
    else
        sd_active = sd
        rd_active = rd
    end

    @info " - Filtering reactions..."
    mask = get_filter_mask(method.filter, sd_active, rd_active)
    splice!(rd_active, findall(mask))
    @info "   - Removed $(count(mask)) filtered reactions from network"

    setup_network!(sd_active, rd_active, method.calculator)
    split_method = method.pars.solve_chunks ? :chunkwise : :complete
    sol = solve_network(method, sd_active, rd_active, Val(split_method), return_integrator)

    if return_integrator
        return sol
    else
        res = ODESolveOutput(method, sol, sd_active, rd_active)
        return res
    end
end

function solve_network(method::StaticODESolve, sd::SpeciesData, rd::RxData, ::Val{:complete}, return_integrator)
    @info " - Removing low-rate reactions"; flush_log()
    apply_low_k_cutoff!(rd, method.calculator, method.pars, method.conditions)

    @info " - Calculating rate constants"; flush_log()
    rates = get_initial_rates(method.conditions, method.calculator)

    @info " - Setting up ReactionSystem"; flush_log()
    @parameters k[1:rd.nr]
    k = collect(k)
    @independent_variables t 
    @species (spec(t))[1:sd.n]
    spec = collect(spec)

    u0 = make_u0(sd, method.pars)
    u0map = Pair.(collect(spec), u0)
    pmap = Pair.(collect(k), rates)

    rs = complete(make_rs(k, spec, t, rd))
    osys = structural_simplify(convert(ODESystem, rs))

    @info " - Formulating ODEProblem"
    @info "   - Sparse? $(method.pars.sparse)"
    @info "   - Analytic Jacobian? $(method.pars.jac)"
    flush_log()
    oprob = ODEProblem(osys, u0map, method.pars.tspan, pmap;
        jac=method.pars.jac, sparse=method.pars.sparse)
    solvecall_kwargs = Dict{Symbol, Any}(
        :progress => method.pars.progress,
        :progress_steps => 10,
        :abstol => method.pars.abstol,
        :reltol => method.pars.reltol,
        :dtmin => eps(method.pars.tspan[end]),
        :maxiters => method.pars.maxiters,
        :saveat => isnothing(method.pars.save_interval) ? eltype(method.pars.tspan)[] : method.pars.save_interval,
        :kwargshandle => KeywordArgError
    )
    if method.pars.ban_negatives
        solvecall_kwargs[:isoutofdomain] = (u,p,t)->any(x->x<0,u)
    end

    @info " - Setting up integrator..."
    integ = init(oprob, method.pars.solver; solvecall_kwargs...)
    if return_integrator
        @info " - Returning integrator early.\n"
        return integ
    end

    adaptive_solve!(integ, method.pars, solvecall_kwargs; print_status=true)

    return integ.sol
end

function solve_network(method::StaticODESolve, sd::SpeciesData, rd::RxData, ::Val{:chunkwise}, return_integrator)
    @info " - Removing low-rate reactions"; flush_log()
    apply_low_k_cutoff!(rd, method.calculator, method.pars, method.conditions)

    @info " - Calculating rate constants"; flush_log()
    rates = get_initial_rates(method.conditions, method.calculator)

    @info " - Setting up ReactionSystem"; flush_log()
    @parameters k[1:rd.nr]
    k = collect(k)
    @independent_variables t 
    @species (spec(t))[1:sd.n]
    spec = collect(spec)

    u0 = make_u0(sd, method.pars)
    u0map = Pair.(collect(spec), u0)
    pmap = Pair.(collect(k), rates)

    rs = complete(make_rs(k, spec, t, rd))
    osys = structural_simplify(convert(ODESystem, rs))
    @info " - Created ReactionSystem"

    @info " - Formulating ODEProblem"
    @info "   - Sparse? $(method.pars.sparse)"
    @info "   - Analytic Jacobian? $(method.pars.jac)"
    flush_log()
    tType = eltype(method.pars.tspan)
    uType = eltype(u0)
    local_tspan = (tType(0.0), tType(method.pars.solve_chunkstep))
    oprob = ODEProblem(osys, u0map, local_tspan, pmap;
        jac=method.pars.jac, sparse=method.pars.sparse)

    # Determine how many solution chunks will be required.
    n_chunks_reqd = Int(method.pars.tspan[2] / method.pars.solve_chunkstep)
    save_interval = isnothing(method.pars.save_interval) ? method.pars.solve_chunkstep : method.pars.save_interval
    saveat_local = collect(0.0:save_interval:method.pars.solve_chunkstep)

    # Allocate final solution arrays.
    size_final = (length(saveat_local)-1)*n_chunks_reqd + 1
    u_final = [zeros(uType, length(u0)) for _ in 1:size_final]
    t_final = zeros(tType, size_final)

    # Initialise integrator.
    solvecall_kwargs = Dict{Symbol, Any}(
        :progress => false,
        :abstol => method.pars.abstol,
        :reltol => method.pars.reltol,
        :dtmin => eps(method.pars.solve_chunkstep),
        :maxiters => method.pars.maxiters,
        :saveat => saveat_local,
        :kwargshandle => KeywordArgError
    )
    if method.pars.ban_negatives
        solvecall_kwargs[:isoutofdomain] = (u,p,t)->any(x->x<0,u)
    end
    @info " - Setting up integrator..."
    integ = init(oprob, method.pars.solver; solvecall_kwargs...)
    if return_integrator
        @info " - Returning integrator early."
        @info "   - Note that integrators for chunkwise solutions require significant work to fully solve outside of their intended solution methods.\n"
        return integ
    end

    # Set up progress bar (if required).
    if method.pars.progress
        pbar_sid = uuid4()
        with_global_logger() do
            @info Progress(pbar_sid, name="Chunkwise ODE")
        end
    end
    @info " - Solving network..."

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
        if nc == n_chunks_reqd-1
            end_idx += 1
            for (i, idx) in enumerate(start_idx:end_idx)
                u_final[idx] .= integ.sol.u[i]
            end
            t_final[start_idx:end_idx] .= integ.sol.t .+ (nc*method.pars.solve_chunkstep)
        # Otherwise, insert all but the final saved timesteps
        else
            for (i, idx) in enumerate(start_idx:end_idx)
                u_final[idx] .= integ.sol.u[i]
            end
            t_final[start_idx:end_idx] .= integ.sol.t[1:end-1] .+ (nc*method.pars.solve_chunkstep)
        end
    end

    if method.pars.progress 
        with_global_logger() do 
            @info Progress(pbar_sid, done=true)
        end 
    end
    @info " - Solved.\n"

    sol = SciMLBase.build_solution(
        oprob,
        method.pars.solver,
        t_final,
        u_final
    )
    return sol
end


"""
    solve_network(method::VariableODESolve, sd::SpeciesData, rd::RxData[, copy_network=true, return_integrator=false])

Solve a network with variable kinetics.

Automatically dispatches to the correct method based
on the value of `method.pars.solve_chunks` and 
`method.conditions.ts_update`, as chunkwise and discrete
solutions require significantly different approaches.

Setting `copy_network=true` generates a `deepcopy` of
the original network in `rd` and `sd` and uses these in
the solution, to avoid side effects from calculators
modifying the original network that is passed in. The
copied (modified) network is returned as part of the 
resulting `ODESolveOutput`.

Setting `return_integrator=true` sets up and returns the
underlying integrator without solving (i.e. at 
`t = method.pars.tspan[1]`), allowing for manual stepping
through the solution. Note that chunkwise solutions implement
many reinitialisations of this integrator, which will have
to be mirrored in the calling script to get the same results.
"""
function solve_network(method::VariableODESolve, sd::SpeciesData, rd::RxData; 
                       copy_network=true, return_integrator=false)
    if copy_network
        sd_active = deepcopy(sd)
        rd_active = deepcopy(rd)
    else
        sd_active = sd
        rd_active = rd
    end

    @info " - Calculating variable condition profiles."; flush_log()
    solve_variable_conditions!(method.conditions, method.pars)

    @info " - Filtering reactions..."
    mask = get_filter_mask(method.filter, sd_active, rd_active)
    splice!(rd_active, findall(mask))
    @info "   - Removed $(count(mask)) filtered reactions from network"

    @info " - Performing calculator-specific network setup."; flush_log()
    setup_network!(sd_active, rd_active, method.calculator)
    split_method = method.pars.solve_chunks ? :chunkwise : :complete
    update_method = method.conditions.discrete_updates ? :discrete : :continuous
    sol = solve_network(method, sd_active, rd_active, Val(split_method), Val(update_method), return_integrator)

    if return_integrator
        return sol
    else
        res = ODESolveOutput(method, sol, sd_active, rd_active)
        return res
    end
end


function solve_network(method::VariableODESolve, sd::SpeciesData, rd::RxData, ::Val{:complete}, ::Val{:continuous}, return_integrator)
    @info " - Removing low-rate reactions"; flush_log()
    apply_low_k_cutoff!(rd, method.calculator, method.pars, method.conditions)

    n_variable_conditions = count(isvariable.(method.conditions.profiles))
    variable_condition_symbols = [sym for sym in method.conditions.symbols if isvariable(method.conditions, sym)]

    @info " - Setting up ReactionSystem"; flush_log()
    @independent_variables t 
    @species (spec(t))[1:sd.n]
    spec = collect(spec)
    @variables (k(t))[1:rd.nr]
    k = collect(k)
    @variables (vc(t))[1:n_variable_conditions]
    vc = collect(vc)

    u0 = make_u0(sd, method.pars)
    u0map = Pair.(collect(spec), u0)
    initial_conditions = get_initial_conditions(method.conditions)
    kmap = Pair.(collect(k), method.calculator(; initial_conditions...))
    vcmap = Pair.(collect(vc), [p.second for p in initial_conditions])
    u0map = vcat(u0map, kmap, vcmap)

    vc_symmap = Dict(sym => num for (sym, num) in zip(variable_condition_symbols, collect(vc)))

    # Form constraint system.
    D = Differential(t)
    direct_profiles = []
    gradient_profiles = []
    gradient_profile_symbols = Symbol[]
    for sym in keys(vc_symmap)
        prof = get_profile(method.conditions, sym)
        if isdirectprofile(prof)
            push!(direct_profiles, vc_symmap[sym] ~ prof.f(t, prof))
        elseif isgradientprofile(prof)
            push!(gradient_profiles, D(vc_symmap[sym]) ~ prof.grad(t, prof))
            push!(gradient_profile_symbols, sym)
        else
            throw(ErrorException("Undefined condition profile type. Something is very wrong..."))
        end
    end
    bound_conditions = vcat(
        [Pair(sym, vc_symmap[sym]) for sym in keys(vc_symmap)], 
        get_static_conditions(method.conditions)
    )
    
    @named rate_sys = ODESystem(
        reduce(vcat, [
            direct_profiles, 
            gradient_profiles, 
            k .~ method.calculator(; bound_conditions...)
        ]), t
    )
    @info "   - Created constraint system for variable conditions."; flush_log()

    rs = make_rs(k, spec, t, rd; variable_k=true)
    @named rs_constrained = extend(rate_sys, rs)
    rs_constrained = complete(rs_constrained)
    @info "   - Merged ReactionSystem with constraints."
    @info "   - Creating ODESystem."
    flush_log()
    osys = structural_simplify(convert(ODESystem, rs_constrained))

    @info " - Formulating ODEProblem"
    @info "   - Sparse? $(method.pars.sparse)"
    @info "   - Analytic Jacobian? $(method.pars.jac)"
    flush_log()
    oprob = ODEProblem(osys, u0map, method.pars.tspan;
        jac=method.pars.jac, sparse=method.pars.sparse,
        warn_initialize_determined=false)
    solvecall_kwargs = Dict{Symbol, Any}(
        :progress => method.pars.progress,
        :progress_steps => 10,
        :abstol => method.pars.abstol,
        :reltol => method.pars.reltol,
        :dtmin => eps(method.pars.tspan[end]),
        :maxiters => method.pars.maxiters,
        :saveat => isnothing(method.pars.save_interval) ? eltype(method.pars.tspan)[] : method.pars.save_interval,
        :tstops => get_tstops(method.conditions),
        :kwargshandle => KeywordArgError
    )
    if method.pars.ban_negatives
        solvecall_kwargs[:isoutofdomain] = (u,p,t)->any(x->x<0,u)
    end

    @info " - Setting up integrator..."
    integ = init(oprob, method.pars.solver; solvecall_kwargs...)
    if return_integrator
        @info " - Returning integrator early.\n"
        return integ
    end

    adaptive_solve!(integ, method.pars, solvecall_kwargs; print_status=true)

    return rebuild_vc_solution(integ.sol, vc_symmap)
end


function solve_network(method::VariableODESolve, sd::SpeciesData, rd::RxData, ::Val{:chunkwise}, ::Val{:continuous}, return_integrator)
    @info " - Removing low-rate reactions"; flush_log()
    apply_low_k_cutoff!(rd, method.calculator, method.pars, method.conditions)

    n_variable_conditions = count(isvariable.(method.conditions.profiles))
    variable_condition_symbols = [sym for sym in method.conditions.symbols if isvariable(method.conditions, sym)]

    @info " - Setting up ReactionSystem"; flush_log()
    @independent_variables t 
    @species (spec(t))[1:sd.n]
    spec = collect(spec)
    @variables (k(t))[1:rd.nr]
    k = collect(k)
    @variables (vc(t))[1:n_variable_conditions]
    vc = collect(vc)
    @parameters chunktime n_chunks

    u0 = make_u0(sd, method.pars)
    u0map = Pair.(collect(spec), u0)
    initial_conditions = get_initial_conditions(method.conditions)
    kmap = Pair.(collect(k), method.calculator(; initial_conditions...))
    vcmap = Pair.(collect(vc), [p.second for p in initial_conditions])
    u0map = vcat(u0map, kmap, vcmap)
    pmap = Pair.([chunktime, n_chunks], [method.pars.solve_chunkstep, 0])

    vc_symmap = Dict(sym => num for (sym, num) in zip(variable_condition_symbols, collect(vc)))

    # Form constraint system.
    D = Differential(t)
    direct_profiles = []
    gradient_profiles = []
    gradient_profile_symbols = Symbol[]
    for sym in keys(vc_symmap)
        prof = get_profile(method.conditions, sym)
        if isdirectprofile(prof)
            push!(direct_profiles, vc_symmap[sym] ~ prof.f(t + (chunktime * n_chunks), prof))
        elseif isgradientprofile(prof)
            push!(gradient_profiles, D(vc_symmap[sym]) ~ prof.grad(t + (chunktime * n_chunks), prof))
            push!(gradient_profile_symbols, sym)
        else
            throw(ErrorException("Undefined condition profile type. Something is very wrong..."))
        end
    end
    bound_conditions = vcat(
        [Pair(sym, vc_symmap[sym]) for sym in keys(vc_symmap)], 
        get_static_conditions(method.conditions)
    )
    
    @named rate_sys = ODESystem(
        reduce(vcat, [
            direct_profiles, 
            gradient_profiles, 
            k .~ method.calculator(; bound_conditions...)
        ]), t
    )
    @info "   - Created constraint system for variable conditions."; flush_log()

    rs = make_rs(k, spec, t, rd, [chunktime, n_chunks])
    @named rs_constrained = extend(rate_sys, rs)
    rs_constrained = complete(rs_constrained)
    @info "   - Merged ReactionSystem with constraints."
    @info "   - Creating ODESystem."
    flush_log()
    osys = structural_simplify(convert(ODESystem, rs_constrained))
    n_chunks_pidx = parameter_index(osys, n_chunks)

    @info " - Formulating ODEProblem"
    @info "   - Sparse? $(method.pars.sparse)"
    @info "   - Analytic Jacobian? $(method.pars.jac)"
    flush_log()
    tType = eltype(method.pars.tspan)
    uType = eltype(u0)
    local_tspan = (0.0, method.pars.solve_chunkstep)
    global_tstops = get_tstops(method.conditions)
    oprob = ODEProblem(osys, u0map, local_tspan, pmap;
        jac=method.pars.jac, sparse=method.pars.sparse,
        warn_initialize_determined=false)

    # Determine how many solution chunks will be required.
    n_chunks_reqd = Int(method.pars.tspan[2] / method.pars.solve_chunkstep)
    save_interval = isnothing(method.pars.save_interval) ? method.pars.solve_chunkstep : method.pars.save_interval
    saveat_local = collect(0.0:save_interval:method.pars.solve_chunkstep)

    # Allocate final solution arrays.
    size_final = (length(saveat_local)-1)*n_chunks_reqd + 1
    u_final = [zeros(uType, length(u0)) for _ in 1:size_final]
    t_final = zeros(tType, size_final)
    vc_final = Dict(sym => zeros(uType, size_final) for sym in gradient_profile_symbols)

    solvecall_kwargs = Dict{Symbol, Any}(
        :progress => false,
        :abstol => method.pars.abstol,
        :reltol => method.pars.reltol,
        :dtmin => eps(method.pars.solve_chunkstep),
        :maxiters => method.pars.maxiters,
        :saveat => saveat_local,
        :kwargshandle => KeywordArgError
    )
    if method.pars.ban_negatives
        solvecall_kwargs[:isoutofdomain] = (u,p,t)->any(x->x<0,u)
    end
    @info " - Setting up integrator..."
    integ = init(oprob, method.pars.solver; solvecall_kwargs...)
    if return_integrator
        @info " - Returning integrator early."
        @info "   - Note that integrators for chunkwise solutions require significant work to fully solve outside of their intended solution methods.\n"
        return integ
    end

    # Set up progress bar (if required).
    if method.pars.progress
        pbar_sid = uuid4()
        with_global_logger() do
            @info Progress(pbar_sid, name="Chunkwise ODE")
        end
    end
    @info " - Solving network..."

    # Loop over the solution chunks needed to generate the full solution.
    for nc in 0:n_chunks_reqd-1
        # Calculate which timestops need to be accounted for in this loop.
        t_start_global = method.pars.solve_chunkstep * nc
        t_end_global = t_start_global + method.pars.solve_chunkstep
        tstops_local = [tg - (nc*method.pars.solve_chunkstep) for tg in global_tstops if (tg >= t_start_global && tg < t_end_global)]

        # Reinitialise the integrator at the current concentrations.
        setindex!(integ.p, nc, n_chunks_pidx)
        chunk_tf = method.pars.solve_chunkstep
        if nc == n_chunks_reqd-1
            # Handle floating point inaccuracies in final timestep.
            last_global_time = tstops_local[end] + nc*method.pars.solve_chunkstep
            if last_global_time == method.pars.tspan[2]
                chunk_tf = tstops_local[end]
            else
                push!(tstops_local, method.pars.solve_chunkstep)
            end
        end
        reinit!(integ, integ.sol.u[end]; tstops=tstops_local, tf=chunk_tf)

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
        if nc == n_chunks_reqd-1
            end_idx += 1
            for (i, idx) in enumerate(start_idx:end_idx)
                u_final[idx] .= integ.sol.u[i][begin:end-length(gradient_profiles)]
            end
            t_final[start_idx:end_idx] .= integ.sol.t .+ (nc*method.pars.solve_chunkstep)
            for (i, sym) in enumerate(gradient_profile_symbols)
                vc_final[sym][start_idx:end_idx] = integ.sol[end-(length(gradient_profiles)-i), :]
            end

        # Otherwise, insert all but the final saved timesteps
        else
            for (i, idx) in enumerate(start_idx:end_idx)
                u_final[idx] .= integ.sol.u[i][begin:end-length(gradient_profiles)]
            end
            t_final[start_idx:end_idx] .= integ.sol.t[1:end-1] .+ (nc*method.pars.solve_chunkstep)
            for (i, sym) in enumerate(gradient_profile_symbols)
                vc_final[sym][start_idx:end_idx] = integ.sol[end-(length(gradient_profiles)-i), 1:end-1]
            end
        end
    end

    if method.pars.progress 
        with_global_logger() do 
            @info Progress(pbar_sid, done=true)
        end 
    end
    @info " - Solved.\n"

    sol = build_vc_solution(
        oprob,
        method.pars.solver,
        t_final,
        u_final,
        vc_final
    )
    return sol
end

function solve_network(method::VariableODESolve, sd::SpeciesData, rd::RxData, ::Val{:complete}, ::Val{:discrete}, return_integrator)
    @info " - Removing low-rate reactions"; flush_log()
    apply_low_k_cutoff!(rd, method.calculator, method.pars, method.conditions)

    @info " - Setting up ReactionSystem"; flush_log()
    @parameters k[1:rd.nr]
    k = collect(k)
    @independent_variables t 
    @species (spec(t))[1:sd.n]
    spec = collect(spec)

    u0 = make_u0(sd, method.pars)
    u0map = Pair.(collect(spec), u0)
    pmap = Pair.(collect(k), method.calculator(; get_initial_conditions(method.conditions)...))

    rs = complete(make_rs(k, spec, t, rd))
    osys = structural_simplify(convert(ODESystem, rs))
    @info " - Created ReactionSystem"

    @info " - Pre-calculating rate constants at discrete time intervals."; flush_log()
    tstops = get_tstops(method.conditions)
    k_precalc = calculate_discrete_rates(method.conditions, method.calculator, rd.nr; uType=eltype(u0))
    k_idx_map = [parameter_index(osys, ki) for ki in k]
    affect! = CompleteRateUpdateAffect(k_precalc, k_idx_map)
    cb = PresetTimeCallback(tstops, affect!; save_positions=(false, false))
    @info " - Created discrete rate constant update callback."

    @info " - Formulating ODEProblem"
    @info "   - Sparse? $(method.pars.sparse)"
    @info "   - Analytic Jacobian? $(method.pars.jac)"
    flush_log()
    oprob = ODEProblem(osys, u0map, method.pars.tspan, pmap;
        jac=method.pars.jac, sparse=method.pars.sparse)
    solvecall_kwargs = Dict{Symbol, Any}(
        :callback => cb,
        :progress => method.pars.progress,
        :progress_steps => 10,
        :abstol => method.pars.abstol,
        :reltol => method.pars.reltol,
        :dtmin => eps(method.pars.tspan[end]),
        :maxiters => method.pars.maxiters,
        :saveat => isnothing(method.pars.save_interval) ? eltype(method.pars.tspan)[] : method.pars.save_interval,
        :tstops => tstops,
        :kwargshandle => KeywordArgError
    )
    if method.pars.ban_negatives
        solvecall_kwargs[:isoutofdomain] = (u,p,t)->any(x->x<0,u)
    end

    @info " - Setting up integrator..."
    integ = init(oprob, method.pars.solver; solvecall_kwargs...)
    if return_integrator
        @info " - Returning integrator early.\n"
        return integ
    end

    adaptive_solve!(integ, method.pars, solvecall_kwargs; print_status=true)

    return build_discrete_rate_solution(integ.sol, k_precalc)
end


function solve_network(method::VariableODESolve, sd::SpeciesData, rd::RxData, ::Val{:chunkwise}, ::Val{:discrete}, return_integrator)
    @info " - Removing low-rate reactions"; flush_log()
    apply_low_k_cutoff!(rd, method.calculator, method.pars, method.conditions)

    @info " - Setting up ReactionSystem"; flush_log()
    @parameters k[1:rd.nr]
    k = collect(k)
    @independent_variables t 
    @species (spec(t))[1:sd.n]
    spec = collect(spec)

    u0 = make_u0(sd, method.pars)
    u0map = Pair.(collect(spec), u0)
    pmap = Pair.(collect(k), method.calculator(; get_initial_conditions(method.conditions)...))

    rs = complete(make_rs(k, spec, t, rd))
    osys = structural_simplify(convert(ODESystem, rs))
    @info " - Created ReactionSystem"

    tType = eltype(method.pars.tspan)
    uType = eltype(u0)
    @info " - Pre-calculating rate constants at discrete time intervals."; flush_log()
    k_precalc = calculate_discrete_rates(method.conditions, method.calculator, rd.nr; uType=uType)
    k_idx_map = [parameter_index(osys, ki) for ki in k]
    condition = ChunkwiseRateUpdateCondition(tType[])
    affect! = ChunkwiseRateUpdateAffect(method.pars.solve_chunkstep, 0, k_precalc, k_idx_map)
    cb = DiscreteCallback(condition, affect!; save_positions=(false, false))
    @info " - Created callback for discrete rate constant updates."

    @info " - Formulating ODEProblem"
    @info "   - Sparse? $(method.pars.sparse)"
    @info "   - Analytic Jacobian? $(method.pars.jac)"
    flush_log()
    local_tspan = (0.0, method.pars.solve_chunkstep)
    global_tstops = get_tstops(method.conditions)
    oprob = ODEProblem(osys, u0map, local_tspan, pmap;
        jac=method.pars.jac, sparse=method.pars.sparse)

    # Determine how many solution chunks will be required.
    n_chunks_reqd = Int(method.pars.tspan[2] / method.pars.solve_chunkstep)
    save_interval = isnothing(method.pars.save_interval) ? method.pars.solve_chunkstep : method.pars.save_interval
    saveat_local = collect(0.0:save_interval:method.pars.solve_chunkstep)

    # Allocate final solution arrays.
    size_final = (length(saveat_local)-1)*n_chunks_reqd + 1
    u_final = [zeros(uType, length(u0)) for _ in 1:size_final]
    t_final = zeros(tType, size_final)

    solvecall_kwargs = Dict{Symbol, Any}(
        :callback => cb,
        :progress => false,
        :abstol => method.pars.abstol,
        :reltol => method.pars.reltol,
        :dtmin => eps(method.pars.solve_chunkstep),
        :maxiters => method.pars.maxiters,
        :saveat => saveat_local,
        :kwargshandle => KeywordArgError
    )
    if method.pars.ban_negatives
        solvecall_kwargs[:isoutofdomain] = (u,p,t)->any(x->x<0,u)
    end
    @info " - Setting up integrator..."
    integ = init(oprob, method.pars.solver; solvecall_kwargs...)
    if return_integrator
        @info " - Returning integrator early."
        @info "   - Note that integrators for chunkwise solutions require significant work to fully solve outside of their intended solution methods.\n"
        return integ
    end

    # Set up progress bar (if required).
    if method.pars.progress
        pbar_sid = uuid4()
        with_global_logger() do
            @info Progress(pbar_sid, name="Chunkwise ODE")
        end
    end
    @info " - Solving network..."

    # Loop over the solution chunks needed to generate the full solution.
    for nc in 0:n_chunks_reqd-1
        # Calculate which timestops need to be accounted for in this loop.
        t_start_global = method.pars.solve_chunkstep * nc
        t_end_global = t_start_global + method.pars.solve_chunkstep
        tstops_local = [tg - (nc*method.pars.solve_chunkstep) for tg in global_tstops if (tg >= t_start_global && tg < t_end_global)]

        # Reinitialise the integrator at the current concentrations.
        condition.tstops_local = tstops_local
        affect!.n_chunks = nc
        chunk_tf = method.pars.solve_chunkstep
        if nc == n_chunks_reqd-1
            # Handle floating point inaccuracies in final timestep.
            last_global_time = tstops_local[end] + nc*method.pars.solve_chunkstep
            if last_global_time == method.pars.tspan[2]
                chunk_tf = tstops_local[end]
            else
                push!(tstops_local, method.pars.solve_chunkstep)
            end
        end
        reinit!(integ, integ.sol.u[end]; tstops=tstops_local, tf=chunk_tf)

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
        if nc == n_chunks_reqd-1
            end_idx += 1
            for (i, idx) in enumerate(start_idx:end_idx)
                u_final[idx] .= integ.sol.u[i]
            end
            t_final[start_idx:end_idx] .= integ.sol.t .+ (nc*method.pars.solve_chunkstep)

        # Otherwise, insert all but the final saved timesteps
        else
            for (i, idx) in enumerate(start_idx:end_idx)
                u_final[idx] .= integ.sol.u[i]
            end
            t_final[start_idx:end_idx] .= integ.sol.t[1:end-1] .+ (nc*method.pars.solve_chunkstep)
        end
    end

    if method.pars.progress 
        with_global_logger() do 
            @info Progress(pbar_sid, done=true)
        end 
    end
    @info " - Solved.\n"

    sol = SciMLBase.build_solution(
        oprob,
        method.pars.solver,
        t_final,
        u_final;
        k = k_precalc,
        retcode = integ.sol.retcode
    )
    return sol
end