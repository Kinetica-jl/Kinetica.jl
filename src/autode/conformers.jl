"""
    autode_conformer_search!(sd::SpeciesData, sid)

Performs a conformer search for the species at `sd.xyz[sid]`, updating its geometry.

Constructs an initial guess of species geometry from its SMILES,
then runs an autodE conformer search with xTB as the energetic
driver. Finds the lowest energy conformer and writes it back
to `sd.xyz[sid]`.

Requires spin multiplicity and charge for the given species to be
cached in `sd.cache[:mult][sid]` and `sd.cache[:charge][sid]` 
respectively, which can be achieved by calling `get_mult!` and 
`get_charge!`. Failure to do so will result in an error.

Also populates `sd.cache` with autodE-derived values for symmetry
number and geometry for later use in TST calculations.
"""
function autode_conformer_search!(sd::SpeciesData, sid)
    if !(sid in keys(sd.cache[:mult])) || !(sid in keys(sd.cache[:charge]))
        throw(KeyError("Missing multiplicity and/or charge in cache for SID $sid."))
    end

    mol = ade.Molecule(smiles=sd.toStr[sid], mult=sd.cache[:mult][sid], charge=sd.cache[:charge][sid])
    if sd.xyz[sid]["N_atoms"] > 2
        @debug "Searching for conformers of species $sid: $(sd.toStr[sid]) (mult = $(sd.cache[:mult][sid]), charge = $(sd.cache[:charge][sid]))"
        mol.find_lowest_energy_conformer()
        n_confs_found = pylen(mol.conformers)
        @debug "$(n_confs_found) conformers found."
    else
        @debug "Species $sid ($(sd.toStr[sid])) too small for conformer search."
        mol.optimise(method=ade.methods.XTB())
    end

    @debug "Writing lowest energy conformer to SpeciesData."
    sd.xyz[sid] = autode_to_frame(mol; info_dict=sd.xyz[sid]["info"])
    sd.xyz[sid]["info"]["energy"] = pyconvert(Float64, mol.energy.to("ev").real)
    
    sd.cache[:symmetry][sid] = pyconvert(Int, mol.symmetry_number)
    if pyconvert(Int, mol.n_atoms) == 1
        sd.cache[:geometry][sid] = 0
    else
        if pyconvert(Bool, mol.is_linear())
            sd.cache[:geometry][sid] = 1
        else
            sd.cache[:geometry][sid] = 2
        end
    end
end

"""
    autode_NCI_conformer_search(sd::SpeciesData, sids[, name="complex"])
    
Performs a search for the lowest energy non-covalent interacting reaction complex defined by the species in `sd` at species IDs `sids`.

Constructs an autodE `NCIComplex` from the provided species
and finds the lowest energy conformation of species. Returns
a new ExtXYZ frame containing the geometry of this conformation.

This frame is additionally tagged with information about the
complex, such as its spin multiplicity, charge and xTB energy.
This information can be found in the resulting `frame`'s "info"
dictionary.

This method should only be performed following calls to 
`autode_conformer_search!`, as it both generates correct species
geometries and also populates `sd.cache` with neccessary
information about the spins and charges of species.

The complexity of this conformer search can be modified by
changing `KineticaASE.ade.Config` before running any calculations,
using the options shown in the autodE documentation
(https://duartegroup.github.io/autodE/examples/nci.html).
"""
function autode_NCI_conformer_search(sd::SpeciesData, sids; name="complex")
    @debug "Searching for conformers of reactive complex."
    mols = [frame_to_autode(sd.xyz[i], mult=sd.cache[:mult][i], chg=sd.cache[:charge][i]) for i in sids]
    sys = ade.NCIComplex(mols..., name=name, do_init_translation=true)

    sys._generate_conformers()
    sys.conformers.optimise(method=ade.methods.XTB())
    sys.conformers.remove_no_energy()
    sys.conformers.prune_diff_graph(graph=sys.graph)
    if pylen(sys.conformers) < 1
        throw(ErrorException("All generated conformers break molecular graph."))
    end
    sys.conformers.prune()
    sys._set_lowest_energy_conformer()
    # sys.find_lowest_energy_conformer(lmethod=ade.methods.XTB())
    n_confs_found = pylen(sys.conformers)
    @debug "$(n_confs_found) conformers found."

    sysframe = autode_to_frame(sys; info_dict=Dict{String, Any}(
        "mult" => pyconvert(Int, sys.mult),
        "chg" => pyconvert(Int, sys.charge),
        "energy" => pyconvert(Float64, sys.energy.to("ev").real),
        "n_species" => length(sids)
    ))
    for f in glob("./*_conf*")
        rm(f)
    end
    return sysframe
end