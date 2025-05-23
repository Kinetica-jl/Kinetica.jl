struct ODESolutionVC{T, N, uType, uType2, DType, tType, rateType, discType, P, A, IType, S,
        AC <: Union{Nothing, Vector{Int}}, R, O, V} <: SciMLBase.AbstractODESolution{T, N, uType}
    u::uType
    u_analytic::uType2
    vcs::Dict{Symbol, Vector{T}}
    errors::DType
    t::tType
    k::rateType
    discretes::discType
    prob::P
    alg::A
    interp::IType
    dense::Bool
    tslocation::Int
    stats::S
    alg_choice::AC
    retcode::ReturnCode.T
    resid::R
    original::O
    saved_subsystem::V
end

function ODESolutionVC{T, N}(u, u_analytic, vcs, errors, t, k, discretes, prob, alg, interp, dense,
    tslocation, stats, alg_choice, retcode, resid, original, saved_subsystem) where {T, N}
    return ODESolutionVC{T, N, typeof(u), typeof(u_analytic), typeof(errors), typeof(t),
                         typeof(k), typeof(discretes), typeof(prob), typeof(alg), typeof(interp),
                         typeof(stats), typeof(alg_choice), typeof(resid), typeof(original), typeof(saved_subsystem)}(
        u, u_analytic, vcs, errors, t, k, discretes, prob, alg, interp,
        dense, tslocation, stats, alg_choice, retcode, resid, original, saved_subsystem)
end


function build_vc_solution(prob::SciMLBase.AbstractODEProblem,
    alg, t, u, vcs;
    dense = false,
    k = nothing,
    alg_choice = nothing,
    interp = SciMLBase.LinearInterpolation(t, u),
    retcode = ReturnCode.Default, destats = missing, stats = nothing,
    kwargs...)

    T = typeof(u[1][1])

    if prob.u0 === nothing
        N = 2
    else
        N = length((size(prob.u0)..., length(u)))
    end

    if typeof(prob.f) <: Tuple
        f = prob.f[1]
    else
        f = prob.f
    end

    if !ismissing(destats)
        msg = "`destats` kwarg has been deprecated in favor of `stats`"
        if stats !== nothing
            msg *= " `stats` kwarg is also provided, ignoring `destats` kwarg."
        else
            stats = destats
        end
        Base.depwarn(msg, :build_solution)
    end

    return ODESolutionVC{T, N}(
        u,
        nothing,
        vcs,
        nothing,
        t, k,
        nothing,
        prob,
        alg,
        interp,
        dense,
        0,
        stats,
        alg_choice,
        retcode,
        nothing, 
        nothing, 
        nothing)
end


function rebuild_vc_solution(sol::ODESolution, vc_symmap::Dict{Symbol, Num})
    sol_params = typeof(sol).parameters
    T = sol_params[1]
    N = sol_params[2]

    vc_dict = Dict{Symbol, Vector{T}}(sym => sol[vc] for (sym, vc) in vc_symmap)

    umat = reduce(vcat, sol.u')[:, begin:end-length(vc_symmap)]
    u = [umat[i, :] for i in axes(umat)[1]]

    return ODESolutionVC{T, N}(
        u,
        nothing,
        vc_dict,
        sol.errors,
        sol.t,
        sol.k,
        sol.discretes,
        sol.prob,
        sol.alg,
        sol.interp,
        sol.dense,
        sol.tslocation,
        sol.stats,
        sol.alg_choice,
        sol.retcode,
        sol.resid,
        sol.original,
        sol.saved_subsystem
    )
end


function build_discrete_rate_solution(sol::ODESolution{T, N}, k::ODESolution{T, N}) where {T, N}
    ODESolution{T, N}(
        sol.u,
        sol.u_analytic,
        sol.errors,
        sol.t,
        k,
        sol.discretes,
        sol.prob,
        sol.alg,
        sol.interp,
        sol.dense,
        sol.tslocation,
        sol.stats,
        sol.alg_choice,
        sol.retcode,
        sol.resid,
        sol.original,
        sol.saved_subsystem)
end

function build_discrete_rate_solution(sol::ODESolution{T, N}, k::DiffEqArray{T, N}) where {T, N}
    ODESolution{T, N}(
        sol.u,
        sol.u_analytic,
        sol.errors,
        sol.t,
        k,
        sol.discretes,
        sol.prob,
        sol.alg,
        sol.interp,
        sol.dense,
        sol.tslocation,
        sol.stats,
        sol.alg_choice,
        sol.retcode,
        sol.resid,
        sol.original,
        sol.saved_subsystem)
end