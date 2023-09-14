module KineticaCore

using Logging
using Dates
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
using PyCall
using PyPlot

const version = VersionNumber(0, 1, 0)

# Global Python package interfaces
const pybel = PyNULL()
const pysys = PyNULL()
const obcr = PyNULL()
function __init__()
    copy!(pybel, pyimport_conda("openbabel.pybel", "openbabel")) # Use pyimport_conda to ensure dependency in place.
    copy!(pysys, pyimport("sys"))
    copy!(obcr, pyimport("obcr"))
end
export pybel, pysys, obcr

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
export NullGradientProfile, LinearGradientProfile
include("conditions/condition_set.jl")
export ConditionSet, isstatic, isvariable
export get_profile, get_tstops, get_t_final
export register_direct_conditions, solve_variable_conditions!

include("openbabel/conversion.jl")
export ingest_xyz_system
include("openbabel/properties.jl")

include("exploration/network.jl")
export SpeciesData, push!, push_unique!
export RxData
include("exploration/cde_utils.jl")
export env_multithread
include("exploration/cde.jl")
export CDE, ingest_cde_run
include("exploration/params.jl")
export DirectExploreParams, IterativeExploreParams
include("exploration/explore_utils.jl")
export import_mechanism, import_mechanism!
include("exploration/methods.jl")
include("exploration/molecule_system.jl")
export system_from_smiles

include("solving/calculator.jl")
export DummyKineticCalculator, PrecalculatedArrheniusCalculator, PrecalculatedLindemannCalculator
export allows_continuous, has_conditions
include("solving/solve_utils.jl")
include("solving/solutions.jl")
include("solving/methods.jl")
export StaticODESolve, VariableODESolve
export solve_network

include("analysis/io.jl")
include("analysis/interpolation.jl")
include("analysis/plotting.jl")

end
