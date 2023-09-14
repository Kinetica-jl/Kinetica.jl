struct ODESolutionVC{T, N, uType, uType2, DType, tType, rateType, P, A, IType, S,
        AC <: Union{Nothing, Vector{Int}}} <: SciMLBase.AbstractODESolution{T, N, uType}
    u::uType
    u_analytic::uType2
    vcs::Dict{Symbol, Vector{T}}
    errors::DType
    t::tType
    k::rateType
    prob::P
    alg::A
    interp::IType
    dense::Bool
    tslocation::Int
    stats::S
    alg_choice::AC
    retcode::ReturnCode.T
end

function ODESolutionVC{T, N}(u, u_analytic, vcs, errors, t, k, prob, alg, interp, dense,
    tslocation, stats, alg_choice, retcode) where {T, N}
    return ODESolutionVC{T, N, typeof(u), typeof(u_analytic), typeof(errors), typeof(t),
        typeof(k), typeof(prob), typeof(alg), typeof(interp),
        typeof(stats),
        typeof(alg_choice)}(u, u_analytic, vcs, errors, t, k, prob, alg, interp,
        dense, tslocation, stats, alg_choice, retcode)
end


function build_vc_solution(prob::SciMLBase.AbstractODEProblem,
    alg, t, u::Vector{Vector{T}}, vcs::Dict{Symbol, Vector{T}};
    dense = false,
    k = nothing,
    alg_choice = nothing,
    interp = SciMLBase.LinearInterpolation(t, u),
    retcode = ReturnCode.Default, destats = missing, stats = nothing,
    kwargs...) where {T}

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

    return ODESolutionVC{T, N}(u,
        nothing,
        vcs,
        nothing,
        t, k,
        prob,
        alg,
        interp,
        dense,
        0,
        stats,
        alg_choice,
        retcode)
end


function rebuild_vc_solution(sol::ODESolution, vc_symbols::Vector{Symbol})
    sol_params = typeof(sol).parameters
    T = sol_params[1]
    N = sol_params[2]

    sym_counter = 0
    vc_dict = Dict{Symbol, Vector{T}}()
    for sym in sol.prob.f.syms
        s = string(sym)
        if occursin("vc", s)
            sym_counter += 1
            vc_dict[vc_symbols[sym_counter]] = sol[sym]
        end
    end
    umat = reduce(vcat, sol.u')[:, begin:end-sym_counter]
    u = [umat[i, :] for i in axes(umat)[1]]

    return ODESolutionVC{T, N}(
        u,
        nothing,
        vc_dict,
        sol.errors,
        sol.t,
        sol.k,
        sol.prob,
        sol.alg,
        sol.interp,
        sol.dense,
        sol.tslocation,
        sol.stats,
        sol.alg_choice,
        sol.retcode
    )
end