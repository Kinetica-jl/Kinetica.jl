"""
Kinetica.jl

UK Ministry of Defence © Crown Owned Copyright 2024/AWE
"""
module Kinetica

using Reexport
using CondaPkg
using Logging
using LoggingExtras
using Dates
using RecipesBase
using Catalyst
@reexport using Catalyst: Graph, savegraph
using OrdinaryDiffEq
using DiffEqCallbacks
using RecursiveArrayTools
using SymbolicIndexingInterface
using LinearAlgebra
using Statistics
using DelimitedFiles
using DataStructures
using ExtXYZ
using ProgressLogging: Progress
using UUIDs: uuid4
using StableHashTraits
using BSON
using Glob
using OrderedCollections
using PythonCall
using CDE_jll

# Global Python package interfaces
const pybel = PythonCall.pynew()
const obcr = PythonCall.pynew()
const pyextxyz = PythonCall.pynew()
const rdChem = PythonCall.pynew()
const rdGeometry = PythonCall.pynew()
const rdLogger = PythonCall.pynew()
const frame_to_rdkit_remap_atoms = PythonCall.pynew()
const np = PythonCall.pynew()
const ase = PythonCall.pynew()
const aseopt = PythonCall.pynew()
const aseneb = PythonCall.pynew()
const aseio = PythonCall.pynew()
const asevib = PythonCall.pynew()
const asethermo = PythonCall.pynew()
const ade = PythonCall.pynew()
const rmsd = PythonCall.pynew()
const rdSmilesParamsWithH = PythonCall.pynew()
function __init__()
    PythonCall.pycopy!(pybel, pyimport("openbabel.pybel"))
    PythonCall.pycopy!(obcr, pyimport("obcr"))
    PythonCall.pycopy!(pyextxyz, pyimport("extxyz"))
    PythonCall.pycopy!(rdChem, pyimport("rdkit.Chem"))
    PythonCall.pycopy!(rdGeometry, pyimport("rdkit.Geometry"))
    PythonCall.pycopy!(rdLogger, pyimport("rdkit.RDLogger"))
    PythonCall.pycopy!(np, pyimport("numpy"))
    PythonCall.pycopy!(ase, pyimport("ase"))
    PythonCall.pycopy!(aseopt, pyimport("ase.optimize"))
    PythonCall.pycopy!(aseneb, pyimport("ase.neb"))
    PythonCall.pycopy!(aseio, pyimport("ase.io"))
    PythonCall.pycopy!(asevib, pyimport("ase.vibrations"))
    PythonCall.pycopy!(asethermo, pyimport("ase.thermochemistry"))
    PythonCall.pycopy!(ade, pyimport("autode"))
    PythonCall.pycopy!(rmsd, pyimport("rmsd"))

    # Set up custom SMILES parser.
    smilesparams = rdChem.SmilesParserParams()
    smilesparams.removeHs = false
    PythonCall.pycopy!(rdSmilesParamsWithH, smilesparams)

    PythonCall.pycopy!(frame_to_rdkit_remap_atoms, pyexec(
        @NamedTuple{f::Py}, 
        """
global pybel; pybel = _pybel
global Chem; Chem = _rdChem

def rdmol_addbonds(pbmol, rdmol):
    for obbond in pybel.ob.OBMolBondIter(pbmol.OBMol):
        a1 = obbond.GetBeginAtom()
        a2 = obbond.GetEndAtom()
        idx1 = a1.GetIdx()
        idx2 = a2.GetIdx()
        rdmol.AddBond(idx1-1, idx2-1, Chem.rdchem.BondType.SINGLE)

    return rdmol.GetMol()

f = rdmol_addbonds""",
        Kinetica,
        (_pybel=pybel, _rdChem=rdChem)
    )[1])

    # Force import of RDKit modules to enable access through rdChem
    pyimport("rdkit.Chem.rdForceFieldHelpers")
    pyimport("rdkit.Chem.rdDistGeom")

    # Disable RdKit logging because it really clogs up the works.
    rdLogger.DisableLog("rdApp.*")

    # Add Conda-installed binaries to PATH
    ENV["PATH"] *= ":"*join(CondaPkg.bindirs(), ":")
end
export pybel, obcr, rdChem

include("constants.jl")
using .Constants

include("logging.jl")
export start_log, end_log, flush_log, flush

include("utils.jl")
export tconvert, create_savepoints
include("pyconvert_utils.jl")

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
export RxData, get_rhash, get_reverse_rhash
export init_network
export format_rxn, print_rxn

include("openbabel/conversion.jl")
export ingest_xyz_system, xyz_to_frame, frame_to_xyz, xyz_file_to_str
export frame_from_smiles, xyz_from_smiles
include("openbabel/properties.jl")
export get_species_stats!

include("rdkit/rdkit.jl")
export frame_to_rdkit, atom_map_smiles, atom_map_frame

include("ase/conversion.jl")
export frame_to_atoms, atoms_to_frame

include("autode/utils.jl")
include("autode/conversion.jl")
include("autode/conformers.jl")

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
include("ase/calculator.jl")
export ASENEBCalculator, calculate_entropy_enthalpy
include("ase/neb.jl")
include("ase/io.jl")
include("ase/optimise.jl")
include("ase/vibrations.jl")
include("ase/builders.jl")
export EMTBuilder, NWChemDFTBuilder, FHIAimsBuilder
include("ase/asethermo_interface.jl")


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
