abstract type AbstractOutputData end

struct ODESolveOutput{sType, kType, vcType} <: AbstractOutputData
    sd::SpeciesData
    rd::RxData
    sol::sType
    sol_k::kType
    sol_vcs::vcType
    pars::ODESimulationParams
    conditions::ConditionSet
end

"""
    ODESolveOutput(solvemethod<:AbstractODESolveMethod, sol<:AbstractODESolution, sd::SpeciesData, rd::RxData)

Data container for output of CRN ODE solutions.

Binds together all data required for analysis of network
ODE solution results. Used as an input for many of the
automated plotting functions in Kinetica.jl.

Also used for results IO, can be completely deconstructed
into a package-independent dictionary tree for saving as
binary JSON (BSON) with `save_output`. Can be reconstructed
from such a BSON file with `load_output`.

Contains fields for:
* CRN `SpeciesData` after simulation (`sd`)
* CRN `RxData` after simulation (`rd`)
* Kinetic simulation solution, usually a DiffEq `ODESolution` o/e (`sol`)
* `DiffEqArray` of precalculated rate constants, if discrete rate update method was used (`sol_k`)
* `Dict` of `DiffEqArray`s for variable conditions, if solved simultaneously with species concentrations (`sol_vcs`)
* `ODESimulationParams` used for kinetic simulation (`pars`)
* `ConditionSet` used for kinetic simulation (`conditions`)
"""
function ODESolveOutput(solvemethod::AbstractODESolveMethod, sol::SciMLBase.AbstractODESolution, sd::SpeciesData, rd::RxData)
    sol_vcs = typeof(sol) <: ODESolutionVC ? Dict(sym => DiffEqArray(val, sol.t) for (sym, val) in sol.vcs) : nothing
    sol_k = typeof(sol.k) <: DiffEqArray ? sol.k : nothing
    return ODESolveOutput(
        sd,
        rd,
        sol,
        sol_k,
        sol_vcs,
        solvemethod.pars,
        solvemethod.conditions
    )
end


"""
    save_output(out::ODESolveOutput, saveto::String)

Saves the output of a CRN ODE solution to BSON file.

Avoids massive file sizes and attempts to maintain forward
compatibility by breaking down large structs into their base
arrays and values, such that they can be reconstructed when
loaded back in.

The resulting BSON file is therefore a dictionary tree, which
can be read in using only Julia's base library and OrderedCollections
if all else fails using `BSON.load` directly.

However, the original `ODESolveOutput` can be mostly reconstructed
by instead calling `load_output()`. Some data is necessarily lost
or converted to a Symbol for reference - mostly data concerning
the internals of `ODESolution`s. 

Additionally, `SurfaceData` structs cannot be saved due to their
many references to Python objects. When saving a surface-based CRN,
the `SpeciesData.surfdata` is therefore not saved, and is set to
`nothing` in the loaded output. 
"""
function save_output(out::ODESolveOutput, saveto::String)
    sol_vcs = !isnothing(out.sol_vcs) ? 
        Dict(sym => val.u for (sym, val) in out.sol_vcs) :
        nothing
    sol_k = typeof(out.sol_k) <: DiffEqArray ? Dict(
                                                    :u => out.sol_k.u,
                                                    :t => out.sol_k.t,
                                               ) : nothing
    condition_profiles = []
    for profile in out.conditions.profiles
        pType = typeof(profile)
        if pType <: AbstractStaticProfile
            pdict = Dict(
                :pType => Symbol(pType),
                :value => profile.value
            )
        elseif pType <: AbstractVariableProfile
            pdict = OrderedDict(fieldnames(pType) .=> getfield.(Ref(profile), fieldnames(pType)))
            if typeof(profile.sol) <: SciMLBase.AbstractODESolution ||
                    typeof(profile.sol) <: RecursiveArrayTools.AbstractDiffEqArray
                pdict[:sol] = Dict(
                    :u => profile.sol.u,
                    :t => profile.sol.t
                )
                pdict[:pType] = Symbol(pType)
            end
            if pType <: AbstractDirectProfile
                pdict[:f] = nothing
            elseif pType <: AbstractGradientProfile
                pdict[:grad] = nothing
            end
        else
            error("Unknown condition profile type for saving.")
        end
        push!(condition_profiles, pdict)
    end

    savedict = Dict(
        :KineticaCoreVersion => pkgversion(Kinetica),
        :sd => Dict(
            :toInt => out.sd.toInt,
            :n => out.sd.n,
            :xyz => out.sd.xyz,
            :level_found => out.sd.level_found
        ),
        :rd => Dict(
            :nr => out.rd.nr,
            :mapped_rxns => out.rd.mapped_rxns,
            :id_reacs => out.rd.id_reacs,
            :id_prods => out.rd.id_prods,
            :stoic_reacs => out.rd.stoic_reacs,
            :stoic_prods => out.rd.stoic_prods,
            :dH => out.rd.dH,
            :rhash => out.rd.rhash,
            :level_found => out.rd.level_found
        ),
        :pars => Dict(
            :tspan => out.pars.tspan,
            :u0 => out.pars.u0,
            :solver => Symbol(typeof(out.pars.solver)),
            :jac => out.pars.jac,
            :sparse => out.pars.sparse,
            :adaptive_tols => out.pars.adaptive_tols,
            :update_tols => out.pars.update_tols,
            :solve_chunks => out.pars.solve_chunks,
            :solve_chunkstep => out.pars.solve_chunkstep,
            :maxiters => out.pars.maxiters,
            :ban_negatives => out.pars.ban_negatives,
            :progress => out.pars.progress,
            :save_interval => out.pars.save_interval,
            :low_k_cutoff => out.pars.low_k_cutoff,
            :allow_short_u0 => out.pars.allow_short_u0
        ),
        :sol => Dict(
            :u => out.sol.u,
            :t => out.sol.t,
            :vcs => sol_vcs,
            :k => sol_k
        ),
        :conditions => Dict(
            :symbols => out.conditions.symbols,
            :profiles => condition_profiles,
            :discrete_updates => out.conditions.discrete_updates,
            :ts_update => out.conditions.ts_update
        )
    )

    bson(saveto, savedict)
end


"""
    load_output(outfile::String)

Loads in the results from a CRN ODE solution generated by `save_output`.

Reconstructs an `ODESolveOutput` from a dictionary tree in
the BSON file `outfile`. Note that some data is lost in the
serialisation process, see the documentation of `save_output`
for details.
"""
function load_output(outfile::String)
    savedict = BSON.load(outfile)
    if savedict[:KineticaCoreVersion] < pkgversion(Kinetica)
        @warn "Loaded network output was made in a previous version of Kinetica.jl, reconstructed output may not be fully compatible."
    elseif savedict[:KineticaCoreVersion] > pkgversion(Kinetica)
        @warn "Loaded network output was made in a newer version of Kinetica.jl, reconstructed output may not be fully compatible."
    end

    sd_iType = typeof(savedict[:sd][:n])
    sd_toStr = Dict{sd_iType, String}(value => key for (key, value) in savedict[:sd][:toInt])
    sd_levels = get(savedict[:sd], :level_found, Dict{sd_iType, Int}(i => 1 for i in 1:savedict[:sd][:n]))
    sd = SpeciesData(
        savedict[:sd][:toInt],
        sd_toStr,
        savedict[:sd][:n],
        savedict[:sd][:xyz],
        nothing,
        sd_levels,
        Dict()
    )
    # For some reason ExtXYZ dicts can get a bit type unstable, species arrays need special handling.
    for i in 1:sd.n
        sd.xyz[i]["arrays"]["species"] = [elem for elem in sd.xyz[i]["arrays"]["species"]]
    end

    rd_iType = typeof(savedict[:rd][:nr])
    rd_mapped_rxns = get(savedict[:rd], :mapped_rxns, String[])
    if length(rd_mapped_rxns) == 0
        @warn "No reaction atom maps found in output."
    end
    rd_levels = get(savedict[:rd], :level_found, ones(Int, savedict[:rd][:nr]))
    rd = RxData(
        savedict[:rd][:nr], String[rxn for rxn in rd_mapped_rxns],
        Vector{rd_iType}[reac for reac in savedict[:rd][:id_reacs]], 
        Vector{rd_iType}[prod for prod in savedict[:rd][:id_prods]],
        Vector{rd_iType}[sreac for sreac in savedict[:rd][:stoic_reacs]], 
        Vector{rd_iType}[sprod for sprod in savedict[:rd][:stoic_prods]],
        savedict[:rd][:dH], 
        Vector{UInt8}[hash for hash in savedict[:rd][:rhash]],
        rd_levels
    )
    if get_rhash(sd, rd, 1) != rd.rhash[1]
        @warn "Reaction hash mismatch, loaded output may not be fully compatible. Consider reassigning reaction hashes."
    end

    pars_fields = fieldnames(ODESimulationParams)
    pars_dict = savedict[:pars]
    for param in keys(pars_dict)
        if !(param in pars_fields)
            pval = pop!(pars_dict, param)
            @warn "Unknown parameter in savedict[:pars] ($param, value: $pval), removing from constructed ODESimulationParams."
        end
    end
    pars = ODESimulationParams(; pars_dict...)

    sol_k = isnothing(savedict[:sol][:k]) ? nothing : 
        DiffEqArray(savedict[:sol][:k][:u], savedict[:sol][:k][:t])
    sol_vcs = isnothing(savedict[:sol][:vcs]) ? nothing :
        Dict(sym => DiffEqArray(val, savedict[:sol][:t]) for (sym, val) in savedict[:sol][:vcs])
    sol = DiffEqArray(savedict[:sol][:u], savedict[:sol][:t])

    profiles = AbstractConditionProfile[]
    for profile_dict in savedict[:conditions][:profiles]
        pType = pop!(profile_dict, :pType)
        if :sol in keys(profile_dict)
            profile_dict[:sol] = DiffEqArray(profile_dict[:sol][:u], profile_dict[:sol][:t])
        end
        if :f in keys(profile_dict)
            profile_dict[:f] = loaded_profile_null_func
        end
        if :grad in keys(profile_dict)
            profile_dict[:grad] = loaded_profile_null_func
        end
        profile = eval(Meta.parse(String(pType)))(values(profile_dict)...)
        push!(profiles, profile)
    end
    conditions = ConditionSet(
        Symbol[sym for sym in savedict[:conditions][:symbols]],
        profiles,
        savedict[:conditions][:discrete_updates],
        savedict[:conditions][:ts_update]
    )

    out = ODESolveOutput(sd, rd, sol, sol_k, sol_vcs, pars, conditions)
    return out
end


function loaded_profile_null_func()
    throw(ErrorException("Condition profile function saving is not supported. Condition function must be manually rebuilt."))
end
loaded_profile_null_func(t) = loaded_profile_null_func()