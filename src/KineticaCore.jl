module KineticaCore

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

include("conditions/variable_temperature.jl")

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
