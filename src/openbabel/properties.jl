"""
    get_species_stats!(sd[, refresh])

Gets statistics about the species in `sd`.

Calculates values for the following useful statistics, which
are placed in the species data cache at `sd.cache` with a key
for each property:

* Number of atoms in species (`:na`)
* Average COM-atom radius of species (`:radii`)
* Molecular weight of species (`:weights`)

Only calculates statistics for species without defined values,
unless `refresh=true` which causes all statistics for every
species to be updated.
"""
function get_species_stats!(sd::SpeciesData{iType}; refresh::Bool=false) where {iType}
    # Create cache dicts if they don't already exist.
    properties = [:radii, :weights]
    property_types = [Float64, Float64]
    for (prop, pType) in zip(properties, property_types)
        if !(prop in keys(sd.cache))
            sd.cache[prop] = Dict{iType, pType}()
        end
    end

    for i in 1:sd.n
        if refresh || !all([i in keys(sd.cache[prop]) for prop in properties])
            path, io = mktemp()
            write_frame(path, sd.xyz[i])
            pbmol = collect(pybel.readfile("xyz", path))[1]

            sd.cache[:weights][i] = pbmol.molwt

            na = sd.xyz[i]["N_atoms"]
            if na == 1
                sd.cache[:radii][i] = pybel.ob.GetVdwRad(pbmol.atoms[1].atomicnum)
            else
                sd.cache[:radii][i] = calc_average_radius(pbmol)
            end
        end
    end
end


"""
    radius = calc_average_radius(pbmol)

Calculates average COM-atom radius of a given Pybel Molecule.

Calculates the radial distance from the center of mass (COM)
to each atom, then takes the mean of these values to obtain
an average molecular radius. Adds on a correction based on the
average van der Walls radius of the atoms to emulate the outer
edge of the molecule.
"""
function calc_average_radius(pbmol)
    na = pbmol.OBMol.NumAtoms()
    masses = zeros(Float64, na)
    coords = zeros(Float64, na, 3)
    atomnums = zeros(Int, na)
    for (j, atom) in enumerate(pbmol.atoms)
        masses[j] = atom.atomicmass
        coords[j, :] = collect(atom.coords)
        atomnums[j] = atom.atomicnum
    end
    molwt = sum(masses)

    com = sum([masses[i] * coords[i, :] for i in 1:na]) / molwt
    r_avg = mean([norm(coords[i, :] - com[:]) for i in 1:na])
    avg_vdw = mean([pybel.ob.GetVdwRad(an) for an in atomnums])
    radius = r_avg + avg_vdw

    return radius
end