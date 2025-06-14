"""
    Graph(sd::SpeciesData, rd::RxData[, graph_attrs, species_attrs, rxn_attrs, edge_attrs, use_smiles=false, remove_inactive_species=true])

Creates a Graphviz graph from the supplied CRN.

Copies the functionality in Catalyst's `ReactionSystem`
graphing utilities (see https://github.com/SciML/Catalyst.jl/blob/master/src/graphs.jl),
but makes them usable on raw Kinetica CRNs. This includes
reworking node names, as Symbolics.jl's 'arrays of symbolic
expressions' are not currently supported in the Catalyst
implementation.

Also extends this functionalty by allowing passing of
graph attributes directly to Graphviz through the 
`graph_attrs` keyword argument. When this is `nothing`,
the default Catalyst graph (drawn using the `dot` layout)
is returned. When given a `Dict{Symbol, String}` of 
graph keywords, these are applied to the generated graph.

Similarly, node attributes for species and reaction nodes
can be passed through `species_attrs` and `rxn_attrs`
respectively. Global edge attributes can be passed through
`edge_attrs`, but individual edges may still set their own
properties.

Additionally allows for plotting with SMILES node labels.
This is disabled by default, as many of the special
characters in SMILES cannot currently be rendered properly.
The exploration level in which species and reactions were
found is also exported as a node attribute.

The resulting `Catalyst.Graph` can be rendered in a notebook,
or saved to file using `Catalyst.savegraph`, both of which
are reexported by Kinetica.
"""
function Catalyst.Graph(sd::SpeciesData, rd::RxData; 
                        graph_attrs::Union{Nothing, Dict{Symbol, String}}=nothing, 
                        species_attrs::Union{Nothing, Dict{Symbol, String}}=nothing,
                        rxn_attrs::Union{Nothing, Dict{Symbol, String}}=nothing,
                        edge_attrs::Union{Nothing, Dict{Symbol, String}}=nothing,
                        use_smiles=false,
                        remove_inactive_species=true)
    gattrs = isnothing(graph_attrs) ? Catalyst.graph_attrs : Catalyst.Attributes(graph_attrs...)
    prog = get(gattrs, :layout, "dot")
    sattrs = isnothing(species_attrs) ? Catalyst.Attributes(:shape => "circle", :color => "#6C9AC3") : Catalyst.Attributes(species_attrs...)
    rattrs = isnothing(rxn_attrs) ? Catalyst.Attributes(:shape => "point", :color => "#E28F41", :width => ".1") : Catalyst.Attributes(rxn_attrs...)
    eattrs = isnothing(edge_attrs) ? Catalyst.edge_attrs : Catalyst.Attributes(edge_attrs...)

    @parameters k[1:rd.nr]
    k = collect(k)
    @independent_variables t 
    @species (spec(t))[1:sd.n]
    spec = collect(spec)
    rs = make_rs(k, spec, t, rd)

    rxs = reactions(rs)
    specs = species(rs)
    spec_names = use_smiles ? [sd.toStr[i] for i in 1:sd.n] : ["S"*subscript(i) for i in 1:sd.n]
    statenodes = []
    states = remove_inactive_species ? unique(reduce(vcat, [rd.id_reacs; rd.id_prods])) : 1:sd.n
    for i in states
        this_sattrs = deepcopy(sattrs)
        this_sattrs[:level] = string(sd.level_found[i])
        push!(statenodes, Catalyst.Node(spec_names[i], this_sattrs))
    end
    transnodes = []
    for (i, r) in enumerate(rxs)
        this_rattrs = deepcopy(rattrs)
        this_rattrs[:level] = string(rd.level_found[i])
        push!(transnodes, Catalyst.Node("R"*subscript(i), this_rattrs))
    end

    stmts = vcat(statenodes, transnodes)
    edges = map(enumerate(rxs)) do (i, r)
        vcat(edgify(zip(r.substrates, r.substoich), spec_names, i, false),
             edgify(zip(r.products, r.prodstoich), spec_names, i, true))
    end
    es = Catalyst.edgifyrates(rxs, specs)
    (!isempty(es)) && push!(edges, es)

    stmts2 = Vector{Catalyst.Statement}()
    append!(stmts2, stmts)
    append!(stmts2, collect(Iterators.flatten(edges)))
    g = Catalyst.Digraph("G", stmts2; prog = prog, graph_attrs = gattrs, 
                         node_attrs = Catalyst.node_attrs, edge_attrs = eattrs)
    return g
end

"""
    edgify(δ, spec_names, i, reverse::Bool)

Rework of Catalyst's `edgify` method to take separate species names.

Symbolics.jl's 'arrays of symbolic expressions' do not
implement the correct fields when passed through Catalyst
`ReactionSystem`s to correctly name graph edges. This
allows for passing an array of species names which match
those used within the main graph, so that edges can be
correctly linked.
"""
function edgify(δ, spec_names, i, reverse::Bool)
    attr = Catalyst.Attributes()
    return map(δ) do p
        val = spec_names[p[1].arguments[2]]
        weight = "$(p[2])"
        attr = Catalyst.Attributes(:label => weight, :labelfontsize => "6")
        return Catalyst.Edge(reverse ? ["R"*subscript(i), "$val"] :
                             ["$val", "R"*subscript(i)], attr)
    end
end

"""
    subscript(i)

Returns the subscripted string of any integer `i`.
"""
subscript(i) = join(Char(0x2080 + d) for d in reverse!(digits(i)))