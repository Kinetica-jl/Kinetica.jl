"""
Kinetica.jl

UK Ministry of Defence © Crown Owned Copyright 2024/AWE
"""
module Kinetica

using Reexport
using Logging
using LoggingExtras
using Dates
using RecipesBase
using Catalyst
@reexport using Catalyst: Graph, savegraph
using OrdinaryDiffEq
using DiffEqCallbacks
using RecursiveArrayTools
using LinearAlgebra
using Statistics
using DelimitedFiles
using DataStructures
using ExtXYZ
using ProgressLogging: Progress
using UUIDs: uuid4
using StableHashTraits
using BSON
using OrderedCollections
using PythonCall
using CDE_jll

const version = VersionNumber(0, 5, 5)

# Global Python package interfaces
const pybel = PythonCall.pynew()
const obcr = PythonCall.pynew()
const pyextxyz = PythonCall.pynew()
const rdChem = PythonCall.pynew()
const rdGeometry = PythonCall.pynew()
const rdLogger = PythonCall.pynew()
function __init__()
    PythonCall.pycopy!(pybel, pyimport("openbabel.pybel"))
    PythonCall.pycopy!(obcr, pyimport("obcr"))
    PythonCall.pycopy!(pyextxyz, pyimport("extxyz"))
    PythonCall.pycopy!(rdChem, pyimport("rdkit.Chem"))
    PythonCall.pycopy!(rdGeometry, pyimport("rdkit.Geometry"))
    PythonCall.pycopy!(rdLogger, pyimport("rdkit.RDLogger"))

    # Force import of RDKit modules to enable access through rdChem
    pyimport("rdkit.Chem.rdForceFieldHelpers")
    pyimport("rdkit.Chem.rdDistGeom")

    # Disable RdKit logging because it really clogs up the works.
    rdLogger.DisableLog("rdApp.*")
end
export pybel, obcr, rdChem

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
export format_rxn, print_rxn

include("openbabel/conversion.jl")
export ingest_xyz_system, xyz_to_frame, frame_to_xyz, xyz_file_to_str
export frame_from_smiles, xyz_from_smiles
include("openbabel/properties.jl")
export get_species_stats!

include("rdkit/rdkit.jl")
export frame_to_rdkit, atom_map_smiles, atom_map_frame

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
include("analysis/graph.jl")

end
