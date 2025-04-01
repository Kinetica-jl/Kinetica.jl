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
function adsorb_frame(frame::Dict{String, Any}, surfdata::SurfaceData, smi::String, heights=nothing)
    surfid = get_surfid(smi)
    siteids = get_surf_siteids(smi)
    if !isnothing(heights) && length(heights) != length(siteids)
        throw(ErrorException("Adsorption heights passed do not match number of adsorbed atoms."))
    end

    surf = surfdata.surfaces[surfid]
    atoms = surf.atoms.copy()
    surf_covradius = mean(
        [pyconvert(Float64, Kinetica.ase.data.covalent_radii[Kinetica.ase.data.atomic_numbers[elem]]) for elem in surf.elements]
    )
    amsmi = atom_map_smiles(frame, smi)
    site_atomids = get_surf_site_atomids(amsmi)

    if length(siteids) == 1
        # Place adsorbate on surface.
        ads_atomid = site_atomids["X$(surfid)_$(siteids[1])"]-1
        ads_atoms = frame_to_atoms(frame)
        ads_atomtype = frame["arrays"]["species"][ads_atomid+1]
        ads_atom_covradius = Kinetica.ase.data.covalent_radii[Kinetica.ase.data.atomic_numbers[ads_atomtype]]
        site_name = surf.sites[siteids[1]]
        height = !isnothing(heights) ? heights[1] : surf_covradius+ads_atom_covradius
        asebuild.add_adsorbate(atoms, ads_atoms, height, site_name, mol_index=ads_atomid)
    else
        # Needs a rigid optimisation that pulls adsorbed atoms to correct sites.
        throw(ErrorException("Adsorption of species at multiple sites is not currently supported."))
    end

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

    surf_atoms = surf.atoms.copy()
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

        # Move gas species to centre of unit cell.
        uc_mid_vector = np.sum(surf_atoms.get_cell()[pyslice(0,2), pyslice(0,2)], axis=0) / 2.0
        gas_positions_centred = gas_positions + np.array([uc_mid_vector[0], uc_mid_vector[1], 0.0])

        # Move bottom of gas species 5 Ang above highest atoms.
        surf_max_z = maximum(combined_atoms.get_positions()[0:pylen(combined_atoms)-1, 2])
        gas_min_z = minimum(gas_positions[0:pylen(gas_atoms)-1, 2])
        zdiff = surf_max_z - gas_min_z
        gas_atoms.set_positions(gas_positions_centred + np.array([0.0, 0.0, zdiff + 5.0 + z_additional]))
        combined_atoms = surf_ads_atoms + gas_atoms

        if has_overlap(combined_atoms)
            # Extend z height of gas species.
            z_additional += 1.0
        elseif !within_unit_cell(combined_atoms)
            # Extend the surface if overlap is detected
            uc_mult += 1
            surf_atoms = surf.atoms.copy().repeat((uc_mult, uc_mult, 1))
        else
            ads_frame = atoms_to_frame(combined_atoms)
            ads_frame["info"]["unit_cell_mult"] = uc_mult
            return ads_frame
        end
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

    surf_atoms = surf.atoms.copy()
    ads1_atoms = frame_to_atoms(ads1)
    ads2_atoms = frame_to_atoms(ads2)
    ads2_pos = ads2_atoms.get_positions()

    attempt = 0
    uc_mult = 1
    while attempt < max_attempts
        attempt += 1

        # Combine adsorbates with the surface
        combined_atoms = surf_atoms + ads1_atoms + ads2_atoms

        if !has_overlap(combined_atoms)
            ads_frame = atoms_to_frame(combined_atoms)
            ads_frame["info"]["unit_cell_mult"] = uc_mult
            return ads_frame
        end

        # Extend the surface if overlap is detected
        uc_mult += 1
        println(uc_mult)
        surf_atoms = surf.atoms.copy().repeat((uc_mult, uc_mult, 1))
        uc_vector = np.sum(surf_atoms.get_cell()[pyslice(0,2), pyslice(0,2)], axis=0)
        ads2_atoms.set_positions(ads2_pos + np.array([uc_vector[0], uc_vector[1], 0.0]))
    end

    throw(ErrorException("Maximum number of attempts reached. Adsorption failed due to overlapping atoms."))
end

"""
    has_overlap(atoms::Py)

Checks if there is any overlap between atoms in the given ASE atoms object.
"""
function has_overlap(atoms::Py)
    positions = pyconvert(Matrix{Float64}, atoms.get_positions())
    for i in axes(positions, 1)[begin:end-1]
        for j in i+1:size(positions, 1)
            if norm(positions[i, :] - positions[j, :]) < 1.0  # Threshold distance for overlap
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