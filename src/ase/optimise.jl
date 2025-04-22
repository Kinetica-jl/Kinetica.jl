"""
    get_mult(smi::String)

Calculates the spin multiplicity of the species represented by SMILES `smi`.

For surfaces-bound species, this is currently undefined and returns
a multiplicity of 1.
"""
get_mult(smi::String) = get_mult(SpeciesStyle(smi), smi)
function get_mult(::GasSpecies, smi::String)
    mol = rdChem.MolFromSmiles(smi)
    n_radical_electrons = 0
    for atom in mol.GetAtoms()
        n_radical_electrons += pyconvert(Int, atom.GetNumRadicalElectrons())
    end
    mult = n_radical_electrons + 1
    return mult
end
get_mult(::SurfaceSpecies, smi::String) = 1

"""
    get_mult(sd::SpeciesData, sid)

Calculates the spin multiplicity of the species in `sd` at ID `sid`.
"""
function get_mult(sd::SpeciesData, sid)
    smi = sd.toStr[sid]
    return get_mult(smi)
end

"""
    get_mult!(sd::SpeciesData, sid)

Caches the spin multiplicity of the species in `sd` at ID `sid` in `sd.cache[:mult]`.
"""
get_mult!(sd::SpeciesData, sid) = sd.cache[:mult][sid] = get_mult(sd, sid)


"""
    get_charge(smi::String)

Calculates the charge of the species represented by SMILES `smi`.

For surfaces-bound species, this is currently undefined and returns
a charge of 0.
"""
get_charge(smi::String) = get_charge(SpeciesStyle(smi), smi)
get_charge(::GasSpecies, smi::String) = pyconvert(Int, pybel.readstring("can", smi).charge)
get_charge(::SurfaceSpecies, smi::String) = 0

"""
    get_charge(sd::SpeciesData, sid)

Calculates the charge of the species in `sd` at ID `sid`.
"""
get_charge(sd::SpeciesData, sid) = get_charge(sd.toStr[sid])

"""
    get_charge!(sd::SpeciesData, sid)

Caches the charge of the species in `sd` at ID `sid` in `sd.cache[:charge]`.
"""
get_charge!(sd::SpeciesData, sid) = sd.cache[:charge][sid] = get_charge(sd, sid)


"""
    get_formal_charges(amsmi::String)
    get_formal_charges(sd::SpeciesData, sid)

Calculates formal charges on each atom of a species.

Can be given an atom-mapped SMILES `amsmi` to calculate from,
or a species in `sd` at ID `sid`, in which case the atom-mapped
SMILES is calculated on-the-fly.

For surfaces-bound species, this is currently undefined and returns
a charge of 0 on each atom.
"""
get_formal_charges(amsmi::String) = get_formal_charges(SpeciesStyle(amsmi), amsmi)
function get_formal_charges(::GasSpecies, amsmi::String)
    mol = pybel.readstring("smi", amsmi)
    formal_charges = [pyconvert(Int, atom.formalcharge) for atom in mol.atoms]
    return formal_charges
end
function get_formal_charges(::SurfaceSpecies, amsmi::String)
    amsmi_replaced = replace(amsmi, r"X\d_\d" => "Fm")
    mol = pybel.readstring("smi", amsmi_replaced)
    n_sites = count(x->x=="Fm", [pyconvert(String, atom.type) for atom in mol.atoms])
    formal_charges = [0 for _ in 1:pylen(mol.atoms)-n_sites]
    return formal_charges
end

function get_formal_charges(sd::SpeciesData, sid)
    amsmi = atom_map_smiles(sd.xyz[sid], sd.toStr[sid])
    return get_formal_charges(amsmi)
end

"""
    get_formal_charges!(sd::SpeciesData, sid)

Caches the formal charges of the species in `sd` at ID `sid` in `sd.cache[:formal_charges]`.
"""
get_formal_charges!(sd::SpeciesData, sid) = sd.cache[:formal_charges][sid] = get_formal_charges(sd, sid)


"""
    get_initial_magmoms(amsmi::String)
    get_initial_magmoms(sd::SpeciesData, sid)

Calculates initial magnetic moments on each atom of a species.

Can be given an atom-mapped SMILES `amsmi` to calculate from,
or a species in `sd` at ID `sid`, in which case the atom-mapped
SMILES is calculated on-the-fly.

For surfaces-bound species, this is currently undefined and returns
a magnetic moment of 0 on each atom.
"""
get_initial_magmoms(amsmi::String) = get_initial_magmoms(SpeciesStyle(amsmi), amsmi)
function get_initial_magmoms(::GasSpecies, amsmi::String)
    mol = rdChem.MolFromSmiles(amsmi, rdSmilesParamsWithH)
    magmoms = zeros(Float64, pyconvert(Int, mol.GetNumAtoms()))
    for atom in mol.GetAtoms()
        magmoms[pyconvert(Int, atom.GetAtomMapNum())] = pyconvert(Float64, atom.GetNumRadicalElectrons())
    end
    return magmoms
end 
function get_initial_magmoms(::SurfaceSpecies, amsmi::String)
    amsmi_replaced = replace(amsmi, r"X\d_\d" => "Fm")
    mol = pybel.readstring("smi", amsmi_replaced)
    n_sites = count(x->x=="Fm", [pyconvert(String, atom.type) for atom in mol.atoms])
    magmoms = [0 for _ in 1:pylen(mol.atoms)-n_sites]
    return magmoms
end

function get_initial_magmoms(sd::SpeciesData, sid)
    amsmi = atom_map_smiles(sd.xyz[sid], sd.toStr[sid])
    return get_initial_magmoms(amsmi)
end

"""
    get_initial_magmoms!(sd::SpeciesData, sid)

Caches the initial magnetic moments of the species in `sd` at ID `sid` in `sd.cache[:initial_magmoms]`.
"""
get_initial_magmoms!(sd::SpeciesData, sid) = sd.cache[:initial_magmoms][sid] = get_initial_magmoms(sd, sid)


"""
    correct_magmoms_for_mult(reac_magmoms::Vector{Float64}, prod_magmoms::Vector{Float64}, mult::Int)

Identifies a set of initial magnetic moments that are consistent with spin multiplicity `mult` across a whole reaction.

As single-reference electronic structure methods cannot 
handle switching of electronic states along a NEB path,
a consistent multiplicity must be used. In cases where
initial magnetic moments are needed, these must be
consistent with the mult for both the reactant and 
product systems.

Attempts to correct existing magmoms by initially flipping lone
radical electrons, e.g. by converting a dissociated system of
two monoradicals to be a spin up-spin down pair with total mult
of 1 to match with a singlet-state reactant. If this cannot be
achieved, resorts to flipping half of an electron pair, e.g. by
converting singlet carbene to triplet carbene.
"""
function correct_magmoms_for_mult!(reac_magmoms::Vector{Float64}, prod_magmoms::Vector{Float64}, mult::Int)
    mdiff(i_magmoms) = (sum(i_magmoms) + 1) - mult

    i_reac_magmoms = [Int(i) for i in reac_magmoms]
    i_prod_magmoms = [Int(i) for i in prod_magmoms]

    # Check for existing match.
    rdiff = mdiff(i_reac_magmoms)
    pdiff = mdiff(i_prod_magmoms)
    if rdiff == 0 && pdiff == 0
        @debug "Initial magnetic moments match reaction spin multiplicity."
        return
    end

    reactive_idxs = [i for i in 1:length(reac_magmoms) if i_reac_magmoms[i] != i_prod_magmoms[i]]
    lone_flippable_reac_idxs = [i for i in reactive_idxs if i_reac_magmoms[i] == 1]
    lone_flippable_prod_idxs = [i for i in reactive_idxs if i_prod_magmoms[i] == 1]
    double_flippable_reac_idxs = [i for i in reactive_idxs if i_reac_magmoms[i] == 2]
    double_flippable_prod_idxs = [i for i in reactive_idxs if i_prod_magmoms[i] == 2]
    if rdiff != 0 && length(lone_flippable_reac_idxs)+length(double_flippable_reac_idxs) == 0
        error("Reactant magmoms cannot be corrected to match reaction multiplicity (no lone radical electrons).")
    elseif pdiff != 0 && length(lone_flippable_prod_idxs)+length(double_flippable_prod_idxs) == 0
        error("Product magmoms cannot be corrected to match reaction multiplicity (no lone radical electrons).")
    end

    while rdiff != 0
        # Prefer lone electron spin flips where possible.
        if length(lone_flippable_reac_idxs) >= abs(rdiff)
            flip_idx = pop!(lone_flippable_reac_idxs)
            @debug "Flipping initial magmom of lone electron in atom $(flip_idx) of reactant."
            i_reac_magmoms[flip_idx] *= -1
            rdiff = mdiff(i_reac_magmoms)
        # Do double flips if diff cannot be resolved with lone flips.
        elseif length(double_flippable_reac_idxs) != 0
            flip_idx = pop!(double_flippable_reac_idxs)
            @debug "Flipping initial magmom of electron pair in atom $(flip_idx) of reactant."
            i_reac_magmoms[flip_idx] = i_reac_magmoms[flip_idx]==0 ? 2 : 0
            rdiff = mdiff(i_reac_magmoms)
        else
            error("Reactant magmoms cannot be corrected to match reaction multiplicity.")
        end
    end

    while pdiff != 0
        # Prefer lone electron spin flips where possible.
        if length(lone_flippable_prod_idxs) >= abs(pdiff)
            flip_idx = pop!(lone_flippable_prod_idxs)
            @debug "Flipping initial magmom of lone electron in atom $(flip_idx) of product."
            i_prod_magmoms[flip_idx] *= -1
            pdiff = mdiff(i_prod_magmoms)
        # Do double flips if diff cannot be resolved with lone flips.
        elseif length(double_flippable_prod_idxs) != 0
            flip_idx = pop!(double_flippable_prod_idxs)
            @debug "Flipping initial magmom of electron pair in atom $(flip_idx) of product."
            i_prod_magmoms[flip_idx] = i_prod_magmoms[flip_idx]==0 ? 2 : 0
            pdiff = mdiff(i_prod_magmoms)
        else
            error("Product magmoms cannot be corrected to match reaction multiplicity.")
        end
    end

    for i in axes(reac_magmoms, 1)
        reac_magmoms[i] = float(i_reac_magmoms[i])
    end
    for i in axes(prod_magmoms, 1)
        prod_magmoms[i] = float(i_prod_magmoms[i])
    end
    return
end


"""
    geomopt!(sd::SpeciesData, i, calc_builder[, calcdir::String="./", optimiser="LBFGSLineSearch", 
             fmax=0.01, maxiters=nothing, check_isomorphic=true, kwargs...])

Runs an ASE-driven geometry optimisation of the species `i` in `sd`.

This is a wrapper around the lower-level `geomopt!(frame, ...)` method
which can dispatch on methods for gas-phase and surface-phase species.
It detects many of the properties needed by the geometry optimisation
from `sd`, and adds extra post-processing of adsorption heights for
the surface-phase case.

`calc_builder` should be a struct with a functor that returns
a correctly constructed ASE calculator for the system at hand.
While ASE can handle many system-specific calculator details
from an `Atoms` object, quantities such as spin multiplicity
and sometimes charge must be input separately. For this reason,
the `calc_builder` functor must take `mult::Int` and `charge::Int`
as its first two arguments. Any other arguments can be passed
via this method's `kwargs`.

This optimisation always uses one of ASE's optimisers out of
"BFGSLineSearch", "fire", "bfgs" or "lbfgs". The maximum force
`fmax` and maximum number of iterations `maxiters` that this
optimiser converges with can also be controlled.

Directly modifies the atomic positions and energy of the
passed in `frame` - for gas-phase species this updates
`sd.xyz[i]`. For surface-phase species, their surface-bound
form gets saved to `sd.cache[:ads_xyz][i]`, while isolated
adsorbates are also updated to `sd.xyz[i]`.

Energies returned are in eV. Returns a boolean for whether the
optimisation converged.
"""
function geomopt!(sd::SpeciesData, i, calc_builder; calcdir::String="./", 
                  optimiser="BFGSLineSearch", fmax=0.01, maxiters=1000, 
                  check_isomorphic=true, kwargs...)
    return geomopt!(SpeciesStyle(sd.toStr[i]), sd, i, calc_builder, calcdir, optimiser, fmax, maxiters,
                    check_isomorphic, kwargs...)
end

function geomopt!(::GasSpecies, sd::SpeciesData, i, calc_builder, calcdir::String, 
                  optimiser, fmax, maxiters, check_isomorphic, kwargs...)
    frame = sd.xyz[i]
    conv = geomopt!(frame, calc_builder; calcdir, mult=sd.cache[:mult][i], chg=sd.cache[:charge][i],
                    formal_charges=sd.cache[:formal_charges][i], initial_magmoms=sd.cache[:initial_magmoms][i],
                    optimiser, fmax, maxiters, check_isomorphic, kwargs...)
    sd.xyz[i] = frame
    return conv
end

function geomopt!(::SurfaceSpecies, sd::SpeciesData, i, calc_builder, calcdir::String, 
                  optimiser, fmax, maxiters, check_isomorphic, kwargs...)
    if isnothing(get(sd.cache, :ads_xyz, nothing))
        throw(ErrorException("Missing cache keys, SpeciesData has not been set up for geometry optimisations."))
    elseif isnothing(get(sd.cache[:ads_xyz], i, nothing))
        throw(ErrorException("Surface species at SID $i is missing an initial conformation, try running `conformer_search!`."))
    end
    frame = sd.cache[:ads_xyz][i]
    copy_frame = deepcopy(frame)
    # Zero these until we have a proper way of defining them for surfaces.
    formal_charges = [0 for _ in 1:frame["N_atoms"]]
    initial_magmoms = [0 for _ in 1:frame["N_atoms"]]
    conv = geomopt!(frame, calc_builder, sd.surfdata; calcdir, mult=sd.cache[:mult][i], chg=sd.cache[:charge][i],
                    formal_charges=formal_charges, initial_magmoms=initial_magmoms,
                    optimiser, fmax, maxiters, check_isomorphic, kwargs...)

    # frame with surface goes to cache - (non)convergence had already been sorted.
    sd.cache[:ads_xyz][i] = frame

    # frame without surface gets information copied over and height tagged, then goes to sd.xyz
    if conv
        atoms = frame_to_atoms(frame)
    else
        atoms = frame_to_atoms(copy_frame)
    end
    _, mol, label = sd.surfdata.finder.predict(atoms)
    if pylen(mol) > 1
        throw(ErrorException("Multiple adsorbates returned from a single-molecule optimisation."))
    end

    adsorbate_frame = atoms_to_frame(mol[0], frame["info"]["energy_ASE"], pyconvert(Vector{Float64}, mol[0].get_moments_of_inertia()))
    adsorbate_frame["info"]["ads_heights"] = [pyconvert(Float64, label[0][i]["height"]) for i in label[0].keys()]
    adsorbate_frame["info"]["adsorbate"] = "true"

    sd.xyz[i] = adsorbate_frame

    return conv
end

"""
    geomopt!(frame::Dict{String, Any}, calc_builder[, calcdir::String="./", mult::Int=1, 
             chg::Int=0, formal_charges=nothing, initial_magmoms=nothing, 
             optimiser="BFGSLineSearch", fmax=0.01, maxiters=1000, check_isomorphic=true, 
             kwargs...]) 
    geomopt!(frame::Dict{String, Any}, calc_builder, surfdata::SurfaceData[, calcdir::String="./",
             mult::Int=1, chg::Int=0, formal_charges=nothing, initial_magmoms=nothing, 
             optimiser="BFGSLineSearch", fmax=0.01, maxiters=1000, check_isomorphic=true, 
             kwargs...]) 

Runs an ASE-driven geometry optimisation of the species in `frame`.

`calc_builder` should be a struct with a functor that returns
a correctly constructed ASE calculator for the system at hand.
While ASE can handle many system-specific calculator details
from an `Atoms` object, quantities such as spin multiplicity
and sometimes charge must be input separately. For this reason,
the `calc_builder` functor must take `mult::Int` and `charge::Int`
as its first two arguments. Any other arguments can be passed
via this method's `kwargs`.

If `surfdata` is provided, this assumes that the system in `frame`
contains molecules on/above a surface and will dispatch a method with
different isomorphism checks.

Some ASE calculators handle charged species at the `Atoms` level.
These require an array of formal charges on each atom to be
given during `Atoms` construction. If `formal_charges` is
provided, this will occur. If not, all formal charges will be 
assumed to be zero. 

IMPORTANTLY, any `formal_charges` array MUST match the atom
ordering of the provided `frame`. This can by calling 
`get_formal_charges` on an atom-mapped SMILES from `atom_map_smiles`.

This optimisation always uses one of ASE's optimisers out of
"BFGSLineSearch", "fire", "bfgs" or "lbfgs". The maximum force
`fmax` and maximum number of iterations `maxiters` that this
optimiser converges with can also be controlled.

Additional graph isomorphism checks will be run through autodE
if `check_isomorphic=true`. If a surface-phase system is
optimised, this checks adosrption sites remain unchanged, and
that each adsorbate maintains its own graph.

Directly modifies the atomic positions and energy of the
passed in `frame`. Energies returned are in eV. Returns a
boolean for whether the optimisation converged.
"""
function geomopt!(frame::Dict{String, Any}, calc_builder; 
                  calcdir::String="./", mult::Int=1, chg::Int=0, 
                  formal_charges=nothing, initial_magmoms=nothing,
                  optimiser="BFGSLineSearch", fmax=0.01, maxiters=1000, 
                  check_isomorphic=true, kwargs...)
    @debug "Starting geometry optimisation."
    atoms = frame_to_atoms(frame, formal_charges, initial_magmoms)
    atoms.calc = calc_builder(calcdir, mult, chg, kwargs...)
    init_energy = pyconvert(Float64, atoms.get_potential_energy())
    init_inertias = pyconvert(Vector{Float64}, atoms.get_moments_of_inertia())

    if optimiser == "BFGSLineSearch"
        opt = aseopt.QuasiNewton(atoms)
    elseif optimiser == "fire"    
        opt = aseopt.FIRE(atoms)
    elseif optimiser == "bfgs"
        opt = aseopt.BFGS(atoms)
    elseif optimiser == "lbfgs"
        opt = aseopt.LBFGS(atoms)
    else
        throw(ArgumentError("Unknown optimiser, must be one of [\"BFGSLineSearch\", \"fire\", \"bfgs\", \"lbfgs\"]"))
    end

    # Optimise with Python exception catching.
    # Also check forces when 10% of the way in to ensure
    # system has not exploded beyond repair.
    conv = false; checkiters = Int(floor(maxiters/10))
    try
        conv = opt.run(fmax=fmax, steps=checkiters)
        conv = pyconvert(Bool, pybuiltins.bool(conv))
        if !conv
            if pyconvert(Float64, opt.get_residual()) > 1e5
                @debug "Optimisation has exploded."
            else
                conv = opt.run(fmax=fmax, steps=maxiters-checkiters)
                conv = pyconvert(Bool, pybuiltins.bool(conv))
            end
        end
    catch err
        conv = false
    end

    if conv && check_isomorphic
        graph_preopt = frame_to_autode(frame; mult=mult, chg=chg).graph
        graph_postopt = frame_to_autode(atoms_to_frame(atoms); mult=mult, chg=chg).graph
        if !autode_is_isomorphic(graph_preopt, graph_postopt)
            conv = false
            @debug "Geometry optimisation breaks molecular graph."
        end
    end

    if conv
        @debug "Geometry optimisation complete."
        frame["arrays"]["pos"] = pyconvert(Matrix, atoms.get_positions().T)
        frame["info"]["energy_ASE"] = pyconvert(Float64, atoms.get_potential_energy())
        frame["arrays"]["inertias"] = pyconvert(Vector{Float64}, atoms.get_moments_of_inertia())
    else
        @debug "Geometry optimisation failed."
        frame["info"]["energy_ASE"] = init_energy
        frame["arrays"]["inertias"] = init_inertias
    end
    return conv
end

function geomopt!(frame::Dict{String, Any}, calc_builder, surfdata::SurfaceData; 
                  calcdir::String="./", mult::Int=1, chg::Int=0, 
                  formal_charges=nothing, initial_magmoms=nothing,
                  optimiser="BFGSLineSearch", fmax=0.01, maxiters=1000, 
                  check_isomorphic=true, kwargs...)
    @debug "Starting geometry optimisation."
    atoms = frame_to_atoms(frame, formal_charges, initial_magmoms)
    _, mols_preopt, labels_preopt = surfdata.finder.predict(atoms)
    atoms.calc = calc_builder(calcdir, mult, chg, kwargs...)
    init_energy = pyconvert(Float64, atoms.get_potential_energy())
    init_inertias = pyconvert(Vector{Float64}, atoms.get_moments_of_inertia())

    if optimiser == "BFGSLineSearch"
        opt = aseopt.QuasiNewton(atoms)
    elseif optimiser == "fire"    
        opt = aseopt.FIRE(atoms)
    elseif optimiser == "bfgs"
        opt = aseopt.BFGS(atoms)
    elseif optimiser == "lbfgs"
        opt = aseopt.LBFGS(atoms)
    else
        throw(ArgumentError("Unknown optimiser, must be one of [\"BFGSLineSearch\", \"fire\", \"bfgs\", \"lbfgs\"]"))
    end

    # Optimise with Python exception catching.
    # Also check forces when 10% of the way in to ensure
    # system has not exploded beyond repair.
    conv = false; checkiters = Int(floor(maxiters/10))
    try
        conv = opt.run(fmax=fmax, steps=checkiters)
        conv = pyconvert(Bool, pybuiltins.bool(conv))
        if !conv
            if pyconvert(Float64, opt.get_residual()) > 1e5
                @debug "Optimisation has exploded."
            else
                conv = opt.run(fmax=fmax, steps=maxiters-checkiters)
                conv = pyconvert(Bool, pybuiltins.bool(conv))
            end
        end
    catch err
        conv = false
    end

    mols_opt = nothing; labels_opt = nothing
    if conv 
        _, mols_opt, labels_opt = surfdata.finder.predict(atoms)
    end

    if conv && check_isomorphic
        # Check that number of molecules is the same.
        # For each molecule, check number of sites is the same.
        # For each site, check its corresponding site in the preopt is the same.
        sites_match = pylen(labels_opt) == pylen(labels_preopt)
        mol_index = 0
        while sites_match && mol_index < pylen(labels_opt)
            sites_match = (pylen(labels_opt[mol_index]) == pylen(labels_preopt[mol_index]))
            if sites_match && pylen(labels_opt[mol_index]) > 0
                sites_match = all([
                    pyconvert(Bool, labels_opt[mol_index][site_index]["site"] == labels_preopt[mol_index][site_index]["site"]) 
                    for site_index in labels_opt[mol_index].keys()
                ])
            end
            mol_index += 1
        end

        if sites_match
            # Use extracted molecules for autodE conformer check.
            for (mol_preopt, mol_opt) in zip(mols_preopt, mols_opt)
                graph_preopt = frame_to_autode(atoms_to_frame(mol_preopt); mult=1, chg=1).graph
                graph_opt = frame_to_autode(atoms_to_frame(mol_opt); mult=1, chg=1).graph
                if !autode_is_isomorphic(graph_preopt, graph_opt)
                    conv = false
                    @debug "Geometry optimisation breaks molecular graph for one or more adsorbates."
                    break
                end
            end
        else
            conv = false
            @debug "Geometry optimisation has moved adsorbate(s) to a different surface site."
        end
    end

    if conv
        @debug "Geometry optimisation complete."
        frame["arrays"]["pos"] = pyconvert(Matrix, atoms.get_positions().T)
        frame["info"]["energy_ASE"] = pyconvert(Float64, atoms.get_potential_energy())
        frame["arrays"]["inertias"] = pyconvert(Vector{Float64}, atoms.get_moments_of_inertia())
        if pylen(mols_opt) == 1
            frame["info"]["ads_heights"] = [pyconvert(Float64, labels_opt[0][i]["height"]) for i in labels_opt[0].keys()]
        end
    else
        @debug "Geometry optimisation failed."
        frame["info"]["energy_ASE"] = init_energy
        frame["arrays"]["inertias"] = init_inertias
    end
    return conv
end


"""
    kabsch_fit!(frame1::Dict{String, Any}, frame2::Dict{String, Any})

Computes the maximum overlap of atoms within `frame1` to `frame2`.

Modifies the atomic positions in `frame1` to be as close as
possible to those in `frame2` through a combination of
translation and rotation. Uses the Kabsch algorithm, as
implemented in the Python package 'rmsd'.
"""
function kabsch_fit!(frame1::Dict{String, Any}, frame2::Dict{String, Any})
    c1 = Py(frame1["arrays"]["pos"]).to_numpy().T
    c2 = Py(frame2["arrays"]["pos"]).to_numpy().T
    frame1["arrays"]["pos"] = pyconvert(Matrix, rmsd.kabsch_fit(c1, c2).T)
    return
end


"""
    get_hydrogen_idxs(amsmi::String)

Returns indices of hydrogen atoms in atom-mapped SMILES `amsmi`.
"""
function get_hydrogen_idxs(amsmi::String)
    at_end = false
    i = 1
    nc = length(amsmi)
    hidxs = [Int[]]
    while !at_end
        if amsmi[i] == '['
            sym = amsmi[i+1]
            idx = parse(Int, string(amsmi[i+3]))
            if sym == 'H'
                push!(hidxs[end], idx)
            end
            i += 5
        elseif amsmi[i] == '.'
            push!(hidxs, Int[])
            i += 1
        else
            i += 1
        end
        if i > nc at_end = true end
    end
    return hidxs
end


"""
    permute_hydrogens!(frame1::Dict{String, Any}, hidxs::Vector{Vector{Int}}, frame2::Dict{String, Any})

Corrects erroneous RDKit atom-mapping on hydrogen atoms.

Due to the implicit handling of hydrogens in RDKit, they can
sometimes receive incorrect atom mapping duing conversion of
a geometry to atom-mapped SMILES. This can cause additional
bond breakage or rotation over the course of a NEB calculation.

This is fixed by permuting every pair of hydrogen atoms in
`frame1` (usually a reactant geometry) and comparing whether
the positional RMSD between it and `frame2` (a product geometry)
decreases. If this is the case, the permutation is made permanent.

This process repeats over all hydrogens in a system until no more
permutations are accepted, indicating the atom mapping with the
least likelihood to cause NEB problems has been found.
"""
function permute_hydrogens!(frame1::Dict{String, Any}, hidxs::Vector{Vector{Int}}, frame2::Dict{String, Any})
    c1 = Py(frame1["arrays"]["pos"]).to_numpy().T
    c2 = Py(frame2["arrays"]["pos"]).to_numpy().T

    if length(reduce(vcat, hidxs)) > 1
        best_pos = c1.copy()
        best_rmsd = pyconvert(Float64, rmsd.kabsch_rmsd(best_pos, c2))
        swapping = true
        while swapping
            has_swapped = false
            for hidxs_mol in hidxs
                if length(hidxs_mol) < 2 continue end
                for i in 1:length(hidxs_mol)-1
                    for j in 2:length(hidxs_mol)
                        swap_pos = best_pos.copy()
                        swap_pos[hidxs_mol[i]-1] = best_pos[hidxs_mol[j]-1]
                        swap_pos[hidxs_mol[j]-1] = best_pos[hidxs_mol[i]-1]
                        swap_rmsd = pyconvert(Float64, rmsd.kabsch_rmsd(swap_pos, c2))
                        if swap_rmsd < best_rmsd
                            @debug "Swapped H$(hidxs_mol[i]) for H$(hidxs_mol[j])"
                            best_pos = swap_pos
                            best_rmsd = swap_rmsd
                            has_swapped = true
                        end
                    end
                end
            end
            if !has_swapped
                swapping = false
            end
        end
        c1 = rmsd.kabsch_fit(best_pos, c2)
    end

    # Any hydrogens with a charge/spin should not isolated and therefore 
    # not swappable, so shouldn't need to be changed here.
    frame1["arrays"]["pos"] = pyconvert(Matrix, c1.T)
    return
end