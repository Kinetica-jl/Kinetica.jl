"""
    autode_get_graph(sd::SpeciesData, sid)

Returns the connectivity graph of a species in `sd` at ID `sid`.

Connectivity graph is returned as a PythonCall `Py` object.
"""
autode_get_graph(sd::SpeciesData, sid) = frame_to_autode(sd.xyz[sid], mult=sd.cache[:mult][sid], chg=sd.cache[:charge][sid]).graph


"""
    autode_is_isomorphic(graph1::Py, graph2::Py)

Returns whether two autodE connectivity graphs are isomorphic.

Graphs can be obtained from `autode_get_graph()`.
"""
autode_is_isomorphic(graph1::Py, graph2::Py) = pyconvert(Bool, ade.mol_graphs.is_isomorphic(graph1, graph2))


"""
    autode_frame_symmetry(frame::Dict{String, Any}[, mult::Int=1, chg::Int=0])

Returns the symmetry number and geometric identifier of a given geometry.

Geometric identifier is 1 if a geometry is linear and 2 otherwise.
This assumes that the geometry is not monoatomic, in which case the
identifier would be 0.
"""
function autode_frame_symmetry(frame::Dict{String, Any}; mult::Int=1, chg::Int=0)
    mol = frame_to_autode(frame; mult=mult, chg=chg)
    sym = pyconvert(Int, mol.symmetry_number)
    if pyconvert(Bool, mol.is_linear())
        geom = 1
    else
        geom = 2
    end
    return sym, geom
end