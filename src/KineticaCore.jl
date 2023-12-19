module KineticaCore

using Logging
using LoggingExtras
using Dates
using RecipesBase
using Catalyst
using DifferentialEquations
using LinearAlgebra
using Statistics
using DelimitedFiles
using DataStructures
using ExtXYZ
using Measurements
using ProgressLogging: Progress
using UUIDs: uuid4
using StableHashTraits
using BSON
using OrderedCollections
using PyCall
using CDE_jll

const version = VersionNumber(0, 2, 0)

# Global Python package interfaces
const pybel = PyNULL()
const pysys = PyNULL()
const obcr = PyNULL()
const pyextxyz = PyNULL()
const rdChem = PyNULL()
const rdLogger = PyNULL()
function __init__()
    copy!(pybel, pyimport_conda("openbabel.pybel", "openbabel")) # Use pyimport_conda to ensure dependency in place.
    copy!(pysys, pyimport("sys"))
    copy!(obcr, pyimport("obcr"))
    copy!(pyextxyz, pyimport("extxyz"))
    copy!(rdChem, pyimport_conda("rdkit.Chem", "rdkit"))
    copy!(rdLogger, pyimport("rdkit.RDLogger"))

    # Disable RdKit logging because it really clogs up the works.
    rdLogger.DisableLog("rdApp.*")
end
export pybel, pysys, obcr, rdChem

include("constants.jl")
using .Constants

include("logging.jl")
export start_log, end_log, flush_log, flush

include("utils.jl")
export tconvert, create_savepoints

include("solving/params.jl")
export ODESimulationParams

include("conditions/abstract_profiles.jl")
export isstatic, isvariable, minimum, maximum
include("conditions/static.jl")
include("conditions/direct_variable.jl")
export NullDirectProfile, LinearDirectProfile
include("conditions/gradient_variable.jl")
export NullGradientProfile, LinearGradientProfile, DoubleRampGradientProfile
include("conditions/condition_set.jl")
export ConditionSet, isstatic, isvariable
export get_profile, get_tstops, get_t_final
export register_direct_conditions, solve_variable_conditions!

include("exploration/network.jl")
export SpeciesData, push!, push_unique!
export RxData
export init_network

include("openbabel/conversion.jl")
export ingest_xyz_system, xyz_to_frame, frame_to_xyz, xyz_file_to_str
include("openbabel/properties.jl")
export get_species_stats!

include("rdkit/rdkit.jl")
export frame_to_rdkit, atom_map_smiles

include("exploration/cde_utils.jl")
export env_multithread
include("exploration/cde.jl")
export CDE, ingest_cde_run
include("exploration/location.jl")
include("exploration/explore_utils.jl")
export import_mechanism, import_mechanism!, import_network
include("exploration/molecule_system.jl")
export system_from_smiles, system_from_mols

include("solving/calculator.jl")
export DummyKineticCalculator, PrecalculatedArrheniusCalculator, PrecalculatedLindemannCalculator
export allows_continuous, has_conditions, setup_network!
include("solving/solve_utils.jl")
include("solving/filters.jl")
export RxFilter, get_filter_mask
include("solving/solutions.jl")
include("solving/methods.jl")
export StaticODESolve, VariableODESolve
export solve_network

include("exploration/methods.jl")
export DirectExplore, IterativeExplore
export explore_network

include("analysis/io.jl")
export ODESolveOutput, save_output, load_output
include("analysis/plotting.jl")

end
