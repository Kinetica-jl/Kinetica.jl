"""
    autode_to_frame(ademol::Py[, info_dict=nothing])

Converts an autodE `Molecule` to an ExtXYZ frame.

Reads atoms from `ademol` and constructs an ExtXYZ frame with
them. Frame defaults to having an empty `info` field, although
this can be populated by passing a `Dict{String, Any}` to the
`info_dict` keyword argument.
"""
function autode_to_frame(ademol::Py; info_dict=nothing)
    na = pylen(ademol.atoms)
    symbols = String[pyconvert(String, ademol.atoms[i].atomic_symbol) for i in 0:na-1]
    coords = reduce(hcat, [
        [
            pyconvert(Float64, ademol.atoms[i].coord.x),
            pyconvert(Float64, ademol.atoms[i].coord.y),
            pyconvert(Float64, ademol.atoms[i].coord.z)
        ] for i in 0:na-1
    ])
    
    info = isnothing(info_dict) ? Dict{String, Any}() : info_dict
    arrays = Dict{String, Any}("species" => symbols, "pos" => coords)
    frame = Dict{String, Any}("N_atoms" => na, "info" => info, "arrays" => arrays)
    return frame
end

"""
    frame_to_autode(frame::Dict{String, Any}[, mult::Int=1, chg::Int=0])

Converts an ExtXYZ frame to an autodE `Molecule`.

Writes a temporary xyz file from `frame` and reads it back in
with autodE. There doesn't seem to be a way around writing to
disk here, as autodE only detects an xyz input when reading a
file ending in '.xyz', which can't be replicated with an in-
memory IO buffer.

Optionally allows specification of spin multiplicity and charge
through the `mult` and `chg` keyword arguments respectively.
"""
function frame_to_autode(frame::Dict{String, Any}; mult::Int=1, chg::Int=0)
    f = joinpath(tempdir(), "frame.xyz")
    write_frame(f, frame)
    mol = ade.Molecule(f, charge=chg, mult=mult)
    rm(f)
    return mol
end
