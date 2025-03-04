"""
    conformer_search!(sd::SpeciesData, sid)

Performs a conformer search for the species at `sd.xyz[sid]`, updating its geometry.

**Gas Phase:** Constructs an initial guess of species geometry from its SMILES,
then runs an autodE conformer search with xTB as the energetic
driver. Finds the lowest energy conformer and writes it back
to `sd.xyz[sid]`.

**Surface Phase:** Uses the starting geometry of an adsorbed species
from `sd.xyz[sid]`, places it on its surface and performs an xTB-based
search of rotamers around the z-axis of the adsorbed atom. Currently
does not support multiply-bound species. Saves the resulting adsorbed
geometry to `sd.cache[:ads_xyz][sid]`.

Requires spin multiplicity and charge for the given species to be
cached in `sd.cache[:mult][sid]` and `sd.cache[:charge][sid]` 
respectively, which can be achieved by calling `get_mult!` and 
`get_charge!`. Failure to do so will result in an error.

Also populates `sd.cache` with autodE-derived values for symmetry
number and geometry for later use in TST calculations.
"""
function conformer_search!(sd::SpeciesData, sid)
    if !(sid in keys(sd.cache[:mult])) || !(sid in keys(sd.cache[:charge]))
        throw(KeyError("Missing multiplicity and/or charge in cache for SID $sid."))
    end
    conformer_search!(SpeciesStyle(sd.toStr[sid]), sd, sid)
    return
end

function conformer_search!(::GasSpecies, sd::SpeciesData, sid)
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

function conformer_search!(::SurfaceSpecies, sd::SpeciesData, sid; n_samples=12)
    print_muter = py_PrintMuter()

    smi = sd.toStr[sid]
    siteids = get_surf_siteids(smi)
    ads_frame = adsorb_frame(sd.xyz[sid], sd.surfdata, smi)

    # Multiply-bound adsorbates can't undergo rotation.
    if length(siteids) == 1
        rotmask = iszero.(ads_frame["arrays"]["tags"]) 
        adsatom_idx = findall(x->x==1, rotmask)[argmin(ads_frame["arrays"]["pos"][3, rotmask])]
        rot_centre = ads_frame["arrays"]["pos"][:, adsatom_idx]
        centre_vecs = [rot_centre for _ in 1:count(rotmask)]

        angles = collect(range(0.0, (2*pi)-(2*pi/n_samples), n_samples))
        v = [[0, 0, 1] for _ in 1:count(rotmask)]
        p = ads_frame["arrays"]["pos"][:, rotmask] .- rot_centre
        pvecs = [p[:, i] for i in axes(p, 2)]

        best_frame = deepcopy(ads_frame)
        best_energy = 0.0
        for a in angles
            rot_frame = deepcopy(ads_frame)
            c = cos(a)
            s = sin(a)
            pvecs_new = c.*pvecs .- cross.(pvecs, s.*v) .+ (dot.(pvecs, v) .* (1.0-c).*v) .+ centre_vecs
            rot_frame["arrays"]["pos"][:, rotmask] = reduce(vcat, pvecs_new')'

            rot_atoms = frame_to_atoms(rot_frame)
            calc = TBLiteBuilder(; method="GFN1-xTB")("./", sd.cache[:mult][sid], sd.cache[:charge][sid])
            rot_atoms.calc = calc
            print_muter.mute()
            rot_energy = pyconvert(Float64, rot_atoms.get_potential_energy())
            print_muter.unmute()
            println(a, " ", rot_energy)

            if rot_energy < best_energy
                best_frame["arrays"]["pos"][:] = rot_frame["arrays"]["pos"][:]
                best_energy = rot_energy
            end
        end

        best_frame["info"]["energy"] = best_energy
        sd.cache[:ads_xyz][sid] = best_frame
    else
        # Just do an xTB energy calculation?
        throw(ErrorException("Multiply-bound adsorbates are not implemented yet."))
    end

    mol = frame_to_autode(sd.xyz[sid]; mult=1, chg=1)
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