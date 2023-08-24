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
using PyPlot

const version = VersionNumber(0, 1, 0)

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

include("solving/solve_utils.jl")
include("solving/methods.jl")

include("analysis/io.jl")
include("analysis/interpolation.jl")
include("analysis/plotting.jl")

end
