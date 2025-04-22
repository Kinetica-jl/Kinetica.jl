"""
    adsorb_frame(sd::SpeciesData, sid::Int[, heights=nothing])
    adsorb_frame(frame::Dict{String, Any}, surfdata::SurfaceData, smi::String[, heights=nothing])

Places the adsorbate `frame` with matching surface SMILES `smi` on its correct surface site.

Reads the desired surface and site from a surface tag in `smi`,
then creates a new frame with the adsorbate on top of the surface
in the correct position.

Defaults to placing an adsorbed atom at a height above the surface
given by the sum of its covalent radius and the mean covalent radius
of the surface elements. Alternatively this can be specified by passing
a vector to `heights` (one element per site).

Currently only suports singly-bound adsorbates, as multiply-bound
adsorbates will need a custom optimisation procedure.
"""
function adsorb_frame(sd::SpeciesData, sid::Int, heights=nothing)
    smi = sd.toStr[sid]
    surfdata = sd.surfdata
    frame = sd.xyz[sid]
    return adsorb_frame(frame, surfdata, smi, heights)
    
end
function adsorb_frame(frame::Dict{String, Any}, surfdata::SurfaceData, smi::String, heights=nothing; max_attempts=10)
    surfid = get_surfid(smi)
    siteids = get_surf_siteids(smi)
    if !isnothing(heights) && length(heights) != length(siteids)
        throw(ErrorException("Adsorption heights passed do not match number of adsorbed atoms."))
    end

    surf = surfdata.surfaces[surfid]
    atoms = pycopy.deepcopy(surf.atoms)
    surf_covradius = mean(
        [pyconvert(Float64, ase.data.covalent_radii[ase.data.atomic_numbers[elem]]) for elem in surf.elements]
    )
    amsmi = atom_map_smiles(frame, smi)
    site_atomids = get_surf_site_atomids(amsmi)

    if length(siteids) == 1
        # Place adsorbate on surface.
        ads_atomid = site_atomids["X$(surfid)_$(siteids[1])"]-1
        ads_atoms = frame_to_atoms(frame)
        ads_atomtype = frame["arrays"]["species"][ads_atomid+1]
        ads_atom_covradius = ase.data.covalent_radii[ase.data.atomic_numbers[ads_atomtype]]
        site_name = surf.sites[siteids[1]]
        height = !isnothing(heights) ? heights[1] : surf_covradius+ads_atom_covradius

        attempts = 0
        uc_mult = 1
        while attempts < max_attempts
            attempts += 1
            surf_na = pylen(atoms)
            asebuild.add_adsorbate(atoms, ads_atoms, height, site_name, mol_index=ads_atomid)

            # Expand surface to check overlaps of periodic adsorbate copies.
            atoms_rep = atoms.repeat((2, 2, 1))
            ads1_idxs = [i+1 for i in surf_na:pylen(atoms)-1]
            ads2_idxs = ads1_idxs .+ pylen(atoms)
            ads3_idxs = ads1_idxs .+ 2*pylen(atoms)
            ads4_idxs = ads1_idxs .+ 3*pylen(atoms)
            if !adsorbed_system_has_overlap(atoms_rep, ads1_idxs, ads2_idxs, 5.0) &&
               !adsorbed_system_has_overlap(atoms_rep, ads1_idxs, ads3_idxs, 5.0) &&
               !adsorbed_system_has_overlap(atoms_rep, ads1_idxs, ads4_idxs, 5.0) &&
               !adsorbed_system_has_overlap(atoms_rep, ads2_idxs, ads3_idxs, 5.0)
                atoms.info["ads_atomid"] = surf_na + ads_atomid
                break
            end

            # If overlaps are detected, repeat the surface and try again.
            uc_mult += 1
            atoms = pycopy.deepcopy(surf.atoms).repeat((uc_mult, uc_mult, 1))
            @debug "Overlaps detected during single-molecule adsorption, extending surface (unit cell multiplier: x$(uc_mult))."
        end

        if attempts >= max_attempts
            throw(ErrorException("Maximum number of attempts reached, adsorption failed."))
        end
    else
        # Needs a rigid optimisation that pulls adsorbed atoms to correct sites.
        throw(ErrorException("Adsorption of species at multiple sites is not currently supported."))
    end

    atoms.info["unit_cell_mult"] = uc_mult
    return atoms_to_frame(atoms)
end

"""
    adsorb_two_frames(sd::SpeciesData, sid1::Int, sid2::Int; max_attempts=10)
    adsorb_two_frames(ads1::Dict{String, Any}, ads2::Dict{String, Any}, surf::Surface, max_attempts::Int=10)

Places two species on a single surface.

Handles the creation of reaction endpoints on surfaces. These can either
be systems with two adsorbates, or a single adsorbate and a gas-phase
species. Double gas-phase species are not allowed, as these should always be
handled by separate elementary reactions. Adsorbates should always include
optimised coordinates and adsorption heights - ensure they have gone through
geomopt! first.

When creating systems with two adsorbates, attempts placement using their
stored coordinates and checks for atom overlaps. If overlaps are detected,
the surface is repeated in the x and y axes and the second adsorbate is
placed in its position within the repeated unit cell. 

When creating mixed adsorbate/gas-phase systems, the adsorbate is placed
using its stored coordinates and the gas-phase species is placed above it 
at a height of 5 Ang. If overlaps are detected, the gas-phase species
is moved up by 1 Ang until no overlaps are detected. If any atoms fall
outside the bounds of the unit cell, e.g. if the gas-phase species is too 
large, the surface is repeated in the x and y axes and the system is recreated.

In either case, if the maximum number of attempts is reached, an error is thrown.
"""
function adsorb_two_frames(sd::SpeciesData, sid1::Int, sid2::Int; max_attempts=10)
    smi1 = sd.toStr[sid1]
    smi2 = sd.toStr[sid2]
    surfid1 = get_surfid(smi1)
    surfid2 = get_surfid(smi2)

    if isnothing(surfid1) && isnothing(surfid2)
        throw(ErrorException("Cannot adsorb two gas-phase species simultaneously."))
    elseif (!isnothing(surfid1) && !isnothing(surfid2)) && (surfid1 != surfid2)
        throw(ErrorException("Cannot adsorb two frames on different surfaces."))
    else
        surfid = isnothing(surfid1) ? surfid2 : surfid1
    end

    surf = sd.surfdata.surfaces[surfid]
    ads1 = sd.xyz[sid1]
    ads2 = sd.xyz[sid2]

    return adsorb_two_frames(XYZStyle(ads1), XYZStyle(ads2), ads1, ads2, surf, max_attempts)
end

function adsorb_two_frames(ads1::Dict{String, Any}, ads2::Dict{String, Any}, surf::Surface, max_attempts::Int=10)
    return adsorb_two_frames(XYZStyle(ads1), XYZStyle(ads2), ads1, ads2, surf, max_attempts)
end
# Handle error for double gas-phase species.
function adsorb_two_frames(::FreeXYZ, ::FreeXYZ, ads1::Dict{String, Any}, ads2::Dict{String, Any}, surf::Surface, max_attempts::Int=10)
    throw(ErrorException("Cannot adsorb two gas-phase species simultaneously."))
end
# Adsorb one molecule, place second above it.
function adsorb_two_frames(::FreeXYZ, ::AdsorbateXYZ, ads1::Dict{String, Any}, ads2::Dict{String, Any}, surf::Surface, max_attempts::Int=10)
    if isnothing(get(ads2["info"], "ads_heights", nothing))
        throw(ErrorException("Adsorption heights not found in adsorbate frame. Please optimise with geomopt! first."))
    end

    surf_atoms = pycopy.deepcopy(surf.atoms)
    gas_atoms = frame_to_atoms(ads1)
    ads_atoms = frame_to_atoms(ads2)

    # Centre gas species at origin.
    gas_positions = gas_atoms.get_positions()
    com = gas_atoms.get_center_of_mass()
    gas_positions -= com

    attempt = 0
    uc_mult = 1
    z_additional = 0.0
    while attempt < max_attempts
        attempt += 1

        # Combine adsorbate with the surface
        surf_ads_atoms = surf_atoms + ads_atoms
        ads1_idxs = [i+1 for i in pylen(surf_atoms):pylen(surf_ads_atoms)-1]

        # Move gas species to centre of unit cell.
        uc_mid_vector = np.sum(surf_atoms.get_cell()[pyslice(0,2), pyslice(0,2)], axis=0) / 2.0
        gas_positions_centred = gas_positions + np.array([uc_mid_vector[0], uc_mid_vector[1], 0.0])

        # Move bottom of gas species 5 Ang above highest atoms.
        surf_max_z = maximum(surf_ads_atoms.get_positions()[0:pylen(surf_ads_atoms)-1, 2])
        gas_min_z = minimum(gas_positions[0:pylen(gas_atoms)-1, 2])
        zdiff = surf_max_z - gas_min_z
        gas_atoms.set_positions(gas_positions_centred + np.array([0.0, 0.0, zdiff + 5.0 + z_additional]))
        combined_atoms = surf_ads_atoms + gas_atoms
        gas1_idxs = [i+1 for i in pylen(surf_ads_atoms):pylen(combined_atoms)-1]

        # Extend z height of gas species if it collides with adsorbate.
        if adsorbed_system_has_overlap(combined_atoms, ads1_idxs, gas1_idxs, 3.0)
            z_additional += 1.0
            @debug "Gas species overlaps with adsorbate in gas-surface system, raising gas molecule by 1 Ang."
            continue
        end

        # Repeat unit cell and check adsorbate and gas species against periodic images.
        atoms_rep = combined_atoms.repeat((2, 2, 1))
        ads2_idxs = ads1_idxs .+ pylen(combined_atoms)
        gas2_idxs = gas1_idxs .+ pylen(combined_atoms)
        ads3_idxs = ads1_idxs .+ 2*pylen(combined_atoms)
        gas3_idxs = gas1_idxs .+ 2*pylen(combined_atoms)
        ads4_idxs = ads1_idxs .+ 3*pylen(combined_atoms)
        gas4_idxs = gas1_idxs .+ 3*pylen(combined_atoms)
        if adsorbed_system_has_overlap(atoms_rep, ads1_idxs, ads2_idxs, 5.0) || # Ads/ads checks
           adsorbed_system_has_overlap(atoms_rep, ads1_idxs, ads3_idxs, 5.0) ||
           adsorbed_system_has_overlap(atoms_rep, ads1_idxs, ads4_idxs, 5.0) ||
           adsorbed_system_has_overlap(atoms_rep, ads2_idxs, ads3_idxs, 5.0) ||
           adsorbed_system_has_overlap(atoms_rep, gas1_idxs, gas2_idxs, 5.0) || # Gas/gas checks
           adsorbed_system_has_overlap(atoms_rep, gas1_idxs, gas3_idxs, 5.0) ||
           adsorbed_system_has_overlap(atoms_rep, gas1_idxs, gas4_idxs, 5.0) ||
           adsorbed_system_has_overlap(atoms_rep, gas2_idxs, gas3_idxs, 5.0) ||
           adsorbed_system_has_overlap(atoms_rep, ads1_idxs, gas2_idxs, 5.0) || # Ads/gas checks
           adsorbed_system_has_overlap(atoms_rep, ads1_idxs, gas3_idxs, 5.0) ||
           adsorbed_system_has_overlap(atoms_rep, ads1_idxs, gas4_idxs, 5.0) ||
           adsorbed_system_has_overlap(atoms_rep, ads2_idxs, gas1_idxs, 5.0) ||
           adsorbed_system_has_overlap(atoms_rep, ads2_idxs, gas3_idxs, 5.0) ||
           adsorbed_system_has_overlap(atoms_rep, ads3_idxs, gas1_idxs, 5.0) ||
           adsorbed_system_has_overlap(atoms_rep, ads3_idxs, gas2_idxs, 5.0) ||
           adsorbed_system_has_overlap(atoms_rep, ads4_idxs, gas1_idxs, 5.0)
            # Extend the surface if overlap is detected
            uc_mult += 1
            surf_atoms = pycopy.deepcopy(surf.atoms).repeat((uc_mult, uc_mult, 1))
            @debug "Overlaps detected during gas-surface adsorption, extending surface (unit cell multiplier: x$(uc_mult))."
            continue
        end

        # If no overlaps are detected, use this configuration.
        ads_frame = atoms_to_frame(combined_atoms)
        ads_frame["info"]["unit_cell_mult"] = uc_mult
        return ads_frame
    end

    throw(ErrorException("Maximum number of attempts reached, adsorption failed."))
end
# Switch order of adsorbates.
function adsorb_two_frames(::AdsorbateXYZ, ::FreeXYZ, ads1::Dict{String, Any}, ads2::Dict{String, Any}, surf::Surface, max_attempts::Int=10)
    return adsorb_two_frames(ads2, ads1, surf, max_attempts)
end
# Adsorb both molecules.
function adsorb_two_frames(::AdsorbateXYZ, ::AdsorbateXYZ, ads1::Dict{String, Any}, ads2::Dict{String, Any}, surf::Surface, max_attempts::Int=10)
    if isnothing(get(ads1["info"], "ads_heights", nothing)) || isnothing(get(ads2["info"], "ads_heights", nothing))
        throw(ErrorException("Adsorption heights not found in adsorbate frames. Please optimise with geomopt! first."))
    end

    surf_atoms = pycopy.deepcopy(surf.atoms)
    ads1_atoms = frame_to_atoms(ads1)
    ads2_atoms = frame_to_atoms(ads2)
    ads2_pos = ads2_atoms.get_positions()

    attempt = 0
    uc_mult = 1
    while attempt < max_attempts
        attempt += 1

        # Combine adsorbates with the surface
        combined_atoms = surf_atoms + ads1_atoms + ads2_atoms
        ads1_idxs = [i+1 for i in pylen(surf_atoms):pylen(surf_atoms)+pylen(ads1_atoms)-1]
        ads2_idxs = [i+1 for i in pylen(surf_atoms)+pylen(ads1_atoms):pylen(combined_atoms)-1]

        # First check that adsorbates in the base cell do not overlap.
        if adsorbed_system_has_overlap(combined_atoms, ads1_idxs, ads2_idxs)
            # Extend the surface if overlap is detected, moving ads2 to match.
            uc_mult += 1
            surf_atoms = pycopy.deepcopy(surf.atoms).repeat((uc_mult, uc_mult, 1))
            uc_vector = np.sum(surf_atoms.get_cell()[pyslice(0,2), pyslice(0,2)], axis=0)/2
            ads2_atoms.set_positions(ads2_pos + np.array([uc_vector[0], uc_vector[1], 0.0]))
            @debug "Overlaps detected during double-adsorbate adsorption, extending surface and moving second adsorbate (unit cell multiplier: x$(uc_mult))."
            continue
        end

        atoms_rep = combined_atoms.repeat((2, 2, 1))
        ads3_idxs = ads1_idxs .+ pylen(combined_atoms)
        ads4_idxs = ads2_idxs .+ pylen(combined_atoms)
        ads5_idxs = ads1_idxs .+ 2*pylen(combined_atoms)
        ads6_idxs = ads2_idxs .+ 2*pylen(combined_atoms)
        ads7_idxs = ads1_idxs .+ 3*pylen(combined_atoms)
        ads8_idxs = ads2_idxs .+ 3*pylen(combined_atoms)

        # Check for overlaps of periodic adsorbate copies.
        if adsorbed_system_has_overlap(atoms_rep, ads1_idxs, ads3_idxs) ||
           adsorbed_system_has_overlap(atoms_rep, ads1_idxs, ads4_idxs) ||
           adsorbed_system_has_overlap(atoms_rep, ads1_idxs, ads5_idxs) ||
           adsorbed_system_has_overlap(atoms_rep, ads1_idxs, ads6_idxs) ||
           adsorbed_system_has_overlap(atoms_rep, ads1_idxs, ads7_idxs) ||
           adsorbed_system_has_overlap(atoms_rep, ads1_idxs, ads8_idxs) ||
           adsorbed_system_has_overlap(atoms_rep, ads2_idxs, ads3_idxs) ||
           adsorbed_system_has_overlap(atoms_rep, ads2_idxs, ads4_idxs) ||
           adsorbed_system_has_overlap(atoms_rep, ads2_idxs, ads5_idxs) ||
           adsorbed_system_has_overlap(atoms_rep, ads2_idxs, ads6_idxs) ||
           adsorbed_system_has_overlap(atoms_rep, ads2_idxs, ads7_idxs) ||
           adsorbed_system_has_overlap(atoms_rep, ads2_idxs, ads8_idxs) ||
           adsorbed_system_has_overlap(atoms_rep, ads3_idxs, ads5_idxs) ||
           adsorbed_system_has_overlap(atoms_rep, ads3_idxs, ads6_idxs) ||
           adsorbed_system_has_overlap(atoms_rep, ads4_idxs, ads5_idxs) ||
           adsorbed_system_has_overlap(atoms_rep, ads4_idxs, ads6_idxs)
            # Extend the surface if overlap is detected, without moving ads2.
            uc_mult += 1
            surf_atoms = pycopy.deepcopy(surf.atoms).repeat((uc_mult, uc_mult, 1))
            @debug "Overlaps detected during double-adsorbate adsorption, extending surface (unit cell multiplier: x$(uc_mult))."
            continue
        end
           
        # If no overlaps are detected, use this configuration.
        ads_frame = atoms_to_frame(combined_atoms)
        ads_frame["info"]["unit_cell_mult"] = uc_mult
        return ads_frame
    end

    throw(ErrorException("Maximum number of attempts reached. Adsorption failed due to overlapping atoms."))
end

"""
    adsorbed_system_has_overlap(atoms::Py, ads1_idxs::Vector{Int}, ads2_idxs::Vector{Int}, min_distance::Float64=2.0)

Checks if there is any overlap between adsorbate atoms in the given ASE atoms object.

Identifies positions of the atoms of each adsorbate and checks
distances between every pair. If any distances are less than
the sum of `min_distance` and the covalent radii of the two atoms,
an overlap is detected.
"""
function adsorbed_system_has_overlap(atoms::Py, ads1_idxs::Vector{Int}, ads2_idxs::Vector{Int}, min_distance::Float64=2.0)
    positions = pyconvert(Matrix{Float64}, atoms.get_positions())
    ads1_pos = positions[ads1_idxs, :]
    ads2_pos = positions[ads2_idxs, :]
    # Get covalent radii of atoms in adsorbates
    ads1_radii = [pyconvert(Float64, Kinetica.ase.data.covalent_radii[Kinetica.ase.data.atomic_numbers[atoms.get_chemical_symbols()[i]]]) for i in ads1_idxs .- 1]
    ads2_radii = [pyconvert(Float64, Kinetica.ase.data.covalent_radii[Kinetica.ase.data.atomic_numbers[atoms.get_chemical_symbols()[i]]]) for i in ads2_idxs .- 1]

    # Check for overlap of atoms across adsorbates, including covalent radii
    for i in axes(ads1_pos, 1)
        for j in axes(ads2_pos, 1)
            if norm(ads1_pos[i, :] - ads2_pos[j, :]) < (ads1_radii[i] + ads2_radii[j] + min_distance)
                return true
            end
        end
    end
    return false
end

"""
    within_unit_cell(frame::Dict{String, Any})
    within_unit_cell(atoms::Py)

Checks if all atoms in the given ASE atoms object are within the unit cell.
"""
function within_unit_cell(frame::Dict{String, Any}, idxs::Vector{Int}=nothing)
    if !isnothing(idxs)
        positions = frame["arrays"]["positions"][1:2, idxs]
    else
        positions = frame["arrays"]["positions"][1:2, :]
    end
    cell = frame["cell"][1:2, 1:2]
    
    # Convert to fractional coordinates
    fractional_positions = inv(cell') * positions
    # Check if all fractional coordinates are within the unit cell
    if any(fractional_positions .< 0) || any(fractional_positions .> 1)
        return false
    end
    return true
end
function within_unit_cell(atoms::Py)
    frame = atoms_to_frame(atoms)
    return within_unit_cell(frame)
end


"""
    center_surface_frame!(frame::Dict{String, Any}, vacuum::Float64=10.0)

Vertically centers atoms in a unit cell with a vacuum above and below.
"""
function center_surface_frame!(frame::Dict{String, Any}, vacuum::Float64=10.0)
    z_min = minimum(frame["arrays"]["pos"][3, :])
    z_max = maximum(frame["arrays"]["pos"][3, :])
    z_diff = z_max - z_min
    frame["cell"][3, 3] = z_diff + 2*vacuum
    frame["arrays"]["pos"][3, :] .+= vacuum - z_min
    return
end


"""
    add_surface_beneath!(frame::Dict{String, Any}, surf::Surface, uc_mult::Int, height=7.5)

Adds a surface beneath the given frame.

The surface is repeated by the given unit cell multiplier `uc_mult` and
placed at the specified height above the frame. The frame is then
centered vertically with a vacuum of 10 Angstroms above and below.
"""
function add_surface_beneath!(frame::Dict{String, Any}, surf::Surface, uc_mult::Int, height=7.5)
    return add_surface_beneath!(XYZStyle(frame), frame, surf, uc_mult, height)
end
function add_surface_beneath!(::AdsorbateXYZ, frame::Dict{String, Any}, surf::Surface, uc_mult::Int, height=7.5)
    throw(ErrorException("Cannot add surface beneath an adsorbate frame. Use adsorb_frame instead."))
end
function add_surface_beneath!(::OnSurfaceXYZ, frame::Dict{String, Any}, surf::Surface, uc_mult::Int, height=7.5)
    throw(ErrorException("Cannot add surface beneath a frame that already contains a surface"))
end
function add_surface_beneath!(::FreeXYZ, frame::Dict{String, Any}, surf::Surface, uc_mult::Int, height=7.5)
    surf_atoms = pycopy.deepcopy(surf.atoms).repeat((uc_mult, uc_mult, 1))
    gas_atoms = frame_to_atoms(frame)

    gas_positions = gas_atoms.get_positions()
    com = gas_atoms.get_center_of_mass()
    gas_positions -= com
    gas_positions[2, :] .+= height
    gas_atoms.set_positions(gas_positions)
    combined_atoms = surf_atoms + gas_atoms
    combined_atoms.center(vacuum=10.0, axis=2)

    combined_frame = atoms_to_frame(combined_atoms)
    frame["N_atoms"] = combined_frame["N_atoms"]
    frame["arrays"] = combined_frame["arrays"]
    frame["cell"] = combined_frame["cell"]
    frame["pbc"] = combined_frame["pbc"]
    frame["info"] = combined_frame["info"]
    frame["info"]["unit_cell_mult"] = uc_mult
    return
end


"""
    scale_surface_to_match!(mod_frame::Dict{String, Any}, ref_frame::Dict{String, Any}, surf::Surface)

Repeats the surface in `mod_frame` to match the size of `ref_frame`.

Determines the unit cell multiplier of the reference frame `ref_frame` using the
base surface in `surf`, then replaces the surface in `mod_frame` with one repeated
by this multiplier. Adsorbates are ignored and left in their original positions.

This is designed to scale surfaces with multipliers greater than 1 - in other words,
`mod_frame` will be extended and not shrunk to match `ref_frame`. If the reference
frame is smaller than the modified frame, the function will throw an error.
"""
function scale_surface_to_match!(mod_frame::Dict{String, Any}, ref_frame::Dict{String, Any}, surf::Surface)
    if ref_frame["N_atoms"] < mod_frame["N_atoms"]
        throw(ErrorException("Reference frame is smaller than modified frame. Cannot scale down."))
    end

    base_atoms = surf.atoms
    
    # Get the unit cell multiplier of the reference frame
    if haskey(ref_frame["info"], "unit_cell_mult")
        ref_uc_mult = ref_frame["info"]["unit_cell_mult"]
    else
        base_cell_xy = pyconvert(Matrix, base_atoms.get_cell())[1:2, 1:2]
        ref_atoms = frame_to_atoms(ref_frame)
        ref_cell_xy = pyconvert(Matrix, ref_atoms.get_cell())[1:2, 1:2]
        ref_uc_mult = round(Int, norm(ref_cell_xy) / norm(base_cell_xy))
    end

    mod_atoms = frame_to_atoms(mod_frame)
    mod_adsatoms = mod_atoms[mod_frame["arrays"]["tags"] .== 0]
    mod_adsatoms.set_cell([0.0, 0.0, 0.0]); mod_adsatoms.set_pbc(false)
    mod_newsurf = pycopy.deepcopy(base_atoms).repeat((ref_uc_mult, ref_uc_mult, 1))
    mod_newatoms = mod_newsurf + mod_adsatoms
    
    mod_newframe = atoms_to_frame(mod_newatoms)
    center_surface_frame!(mod_newframe)
    mod_frame["N_atoms"] = mod_newframe["N_atoms"]
    mod_frame["arrays"] = mod_newframe["arrays"]
    mod_frame["cell"] = mod_newframe["cell"]
    mod_frame["info"] = mod_newframe["info"]
    mod_frame["info"]["unit_cell_mult"] = ref_uc_mult
    return    
end