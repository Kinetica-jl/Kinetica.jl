abstract type AbstractExploreParams end

"""
Container for parameters used in direct CRN exploration.

Should usually be instantiated using one of its constructors,
as CDE should either be instantiated manually or its parameters
should be passed through here, not both.

Contains fields for:
* Top level of CRN exploration directory (`rdir_head`)
* SMILES string(s) of main breakdown reactant(s) being studied (`reac_smiles`)
* Maximum number of iterations to perform (`maxiters`)
* Number of iterations with no change in reactions to consider as converged (`reacthresh`)
* CDE template directory (`cde_tdir`)
* CDE initial breakdown xyz file (`cde_init_xyz`)
* CDE radius to explore in reaction space per mechanism (`cde_radius`)
* CDE number of mechanisms to generate per iteration (`cde_nrxn`)
* Number of parallel CDE runs to perform, must be ≥ `cde_nrxn` (`cde_threads`)
* Number of threads to run xTB calculations with (`xtb_threads`)
* `CDE` instance (`cde`)
"""
struct DirectExploreParams <: AbstractExploreParams
    rdir_head::String
    reac_smiles::Union{String, Vector{String}}
    maxiters::Integer
    reacthresh::Integer

    # CDE params
    cde_tdir::String
    cde_init_xyz::String
    cde_radius::Integer
    cde_nrxn::Integer
    cde_threads::Integer
    xtb_threads::Integer

    # Actual CDE runner
    cde::CDE
end

"""
    epars = DirectExploreParams(;
        rdir_head, reac_smiles, maxiters, reacthresh,
        cde_tdir, cde_init_xyz[, cde_radius, cde_nrxn,
        cde_threads, xtb_threads])

Outer constructor method for `DirectExploreParams`.

Can be used for simple configuration of direct CRN
exploration. Does not require setting up CDE interface 
manually.
"""
function DirectExploreParams(; 
        rdir_head::String,
        reac_smiles::Union{String, Vector{String}},
        maxiters::Integer,
        reacthresh::Integer,
        cde_tdir::String,
        cde_init_xyz::String,
        cde_radius::Integer = 50,
        cde_nrxn::Integer = 1,
        cde_threads::Integer = 1,
        xtb_threads::Integer = 4)
    env = env_multithread(xtb_threads)
    cde = CDE(
        rdir = rdir_head,
        tdir = cde_tdir,
        init_xyz = cde_init_xyz,
        env = env,
        radius = cde_radius,
        nrxn = cde_nrxn,
        parallel_runs = cde_threads,
        parallel_exes = cde_threads
    )
    return DirectExploreParams(
        rdir_head = rdir_head,
        reac_smiles = reac_smiles,
        maxiters = maxiters,
        reacthresh = reacthresh,
        cde_tdir = cde_tdir, 
        cde_init_xyz = cde_init_xyz,
        cde_radius = cde_radius,
        cde_nrxn = cde_nrxn,
        cde_threads = cde_threads,
        xtb_threads = xtb_threads,
        cde = cde
    )
end

"""
    epars = DirectExploreParams(;
        rdir_head, reac_smiles, maxiters, reacthresh, cde)

Outer constructor method for `DirectExploreParams`.

Requires external manual configuration of CDE interface,
which can be more difficult but offers additional customisation.
"""
function DirectExploreParams(;
        rdir_head::String,
        reac_smiles::Union{String, Vector{String}},
        maxiters::Integer,
        reacthresh::Integer,
        cde::CDE)
    return DirectExploreParams(
        rdir_head = rdir_head,
        reac_smiles = reac_smiles,
        maxiters = maxiters,
        reacthresh = reacthresh,
        cde_tdir = cde.tdir, 
        cde_init_xyz = cde.init_xyz,
        cde_radius = cde.radius,
        cde_nrxn = cde.nrxn,
        cde_threads = cde.parallel_exes,
        xtb_threads = cde.env["OMP_NUM_THREADS"],
        cde = cde
    )
end


"""
Container for parameters used in iterative kinetics-based CRN exploration.

Should usually be instantiated using one of its constructors,
as CDE should either be instantiated manually or its parameters
should be passed through here, not both.

Contains fields for:
* Top level of CRN exploration directory (`rdir_head`)
* SMILES string(s) of main breakdown reactant(s) being studied (`reac_smiles`)
* Maximum number of iterations to perform (`maxiters`)
* Number of iterations with no change in reactions to consider as converged (`reacthresh`)
* Concentration above which species will be selected as seeds each level (`seed_conc`)
* Blacklist of species to avoid doing independent supspace explorations on (`independent_blacklist`)
* CDE template directory (`cde_tdir`)
* CDE initial breakdown xyz file (`cde_init_xyz`)
* CDE radius to explore in reaction space per mechanism (`cde_radius`)
* CDE number of mechanisms to generate per iteration (`cde_nrxn`)
* Number of parallel CDE runs to perform, must be ≥ `cde_nrxn` (`cde_threads`)
* Number of threads to run xTB calculations with (`xtb_threads`)
* `CDE` instance (`cde`)
"""
Base.@kwdef struct IterativeExploreParams{uType} <: AbstractExploreParams
    rdir_head::String
    reac_smiles::Union{String, Vector{String}}
    maxiters::Integer
    reacthresh::Integer
    seed_conc::uType
    independent_blacklist::Vector{String}

    # CDE params
    cde_tdir::String
    cde_init_xyz::String
    cde_radius::Integer
    cde_nrxn::Integer
    cde_threads::Integer
    xtb_threads::Integer

    # Actual CDE runner
    cde::CDE
end

"""
    epars = IterativeExploreParams(;
        rdir_head, reac_smiles, maxiters, reacthresh,
        cde_tdir, cde_init_xyz[, seed_conc, 
        independent_blacklist, cde_radius, cde_nrxn,
        cde_threads, xtb_threads])

Outer constructor method for `IterativeExploreParams`.

Can be used for simple configuration of iterative CRN
exploration. Does not require setting up CDE interface 
manually.
"""
function IterativeExploreParams(; 
        rdir_head::String,
        reac_smiles::Union{String, Vector{String}},
        maxiters::Integer,
        reacthresh::Integer,
        cde_tdir::String,
        cde_init_xyz::String,
        seed_conc::uType = 0.05,
        independent_blacklist::Vector{String} = [],
        cde_radius::Integer = 50,
        cde_nrxn::Integer = 1,
        cde_threads::Integer = 1,
        xtb_threads::Integer = 4) where {uType <: AbstractFloat}
    env = env_multithread(xtb_threads)
    cde = CDE(
        rdir = rdir_head,
        tdir = cde_tdir,
        init_xyz = cde_init_xyz,
        env = env,
        radius = cde_radius,
        nrxn = cde_nrxn,
        parallel_runs = cde_threads,
        parallel_exes = cde_threads
    )
    return IterativeExploreParams(
        rdir_head = rdir_head,
        reac_smiles = reac_smiles,
        maxiters = maxiters,
        reacthresh = reacthresh,
        seed_conc = seed_conc,
        independent_blacklist = independent_blacklist,
        cde_tdir = cde_tdir, 
        cde_init_xyz = cde_init_xyz,
        cde_radius = cde_radius,
        cde_nrxn = cde_nrxn,
        cde_threads = cde_threads,
        xtb_threads = xtb_threads,
        cde = cde
    )
end

"""
    epars = IterativeExploreParams(;
        rdir_head, reac_smiles, maxiters, reacthresh, 
        cde[, seed_conc, independent_blacklist])

Outer constructor method for `IterativeExploreParams`.

Requires external manual configuration of CDE interface,
which can be more difficult but offers additional customisation.
"""
function IterativeExploreParams(;
        rdir_head::String,
        reac_smiles::Union{String, Vector{String}},
        maxiters::Integer,
        reacthresh::Integer,
        cde::CDE,
        seed_conc::uType = 0.05,
        independent_blacklist::Vector{String} = [])
    return IterativeExploreParams(
        rdir_head = rdir_head,
        reac_smiles = reac_smiles,
        maxiters = maxiters,
        reacthresh = reacthresh,
        seed_conc = seed_conc,
        independent_blacklist = independent_blacklist,
        cde_tdir = cde.tdir, 
        cde_init_xyz = cde.init_xyz,
        cde_radius = cde.radius,
        cde_nrxn = cde.nrxn,
        cde_threads = cde.parallel_exes,
        xtb_threads = cde.env["OMP_NUM_THREADS"],
        cde = cde
    )
end