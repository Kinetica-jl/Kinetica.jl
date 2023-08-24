module KineticaCore

using Logging
using Dates
using Catalyst
using DifferentialEquations
using LinearAlgebra
using DelimitedFiles
using DataStructures
using ExtXYZ
using Measurements
using ProgressLogging: Progress
using UUIDs: uuid4
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

include("conditions/abstract_profiles.jl")
include("conditions/variable_temperature.jl")
export NullTprofile, LinearTprofile
export SimpleDoubleRampTprofile, SmoothDoubleRampTprofile

include("exploration/network.jl")
include("exploration/explore_utils.jl")
include("exploration/methods.jl")
include("exploration/molecule_system.jl")
export system_from_smiles

include("solving/solve_utils.jl")
include("solving/methods.jl")

include("analysis/io.jl")
include("analysis/interpolation.jl")
include("analysis/plotting.jl")

end
