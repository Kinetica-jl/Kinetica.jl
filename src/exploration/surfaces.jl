mutable struct Surface
    name::String
    atoms::Py
    elements::Set{String}
    sites::Dict{Int, String}
    siteids::Dict{String, Int}
    sitecoords::Dict{String, Int}
    sitebounds::Dict{String, <:Any}
    cutoff_mult::Float64
end

"""
    Surface(name::String, atoms::Py[, cutoff_mult=1.0, sitecoords::Dict{String, Int}=Dict(), sitebounds::Dict{String, <:Any}=Dict()])

Constructs a `Surface` directly from an ASE surface slab.

Can take any `name`, this is just a user-facing identifier for
the surface.

The surface slab supplied in `atoms` should be created with
one of the surface builders in `ase.build` (e.g. `fcc111()`),
such that they are compatible with ASE functions such as 
`ase.build.add_adsorbate()`. 

If `sitebounds` is passed, it should be a `Dict` for string-form 
absorption site names bound to instances of `asesurfacefinder.SampleBounds`.
Details can be found in ASESurfaceFinder's docs. If bounds for any
site are not provided, they will be set to ASESurfaceFinder's defaults
when passing this `Surface` to `SurfaceData`.
"""
function Surface(name::String, atoms::Py; cutoff_mult=1.0, 
                 sitecoords::Dict{String, Int}=Dict{String, Int}(), sitebounds::Dict{String, <:Any}=Dict{String, Any}())
    # Check Atoms object is usable.
    if !pyconvert(Bool, (pytype(atoms) == ase.Atoms))
        error("Cannot construct Surface '$(name)': provided `atoms` are not of type `ase.Atoms`.")
    end
    if !("adsorbate_info" in atoms.info)
        error("Cannot construct Surface '$(name)': no adsorbate info.")
    end

    # Extract site names and bind to integer IDs.
    # Canonicalise order of site ids by sorting names alphabetically.
    site_names = sort(pyconvert(Vector{String}, atoms.info["adsorbate_info"]["sites"].keys()))
    sites = Dict{Int, String}(i => s for (i, s) in enumerate(site_names))
    siteids = Dict{String, Int}(s => i for (i, s) in enumerate(site_names))
    for site in keys(sitecoords)
        if !(site in keys(siteids))
            error("Site '$(site)' provided to sitecoords does not exist on Surface '$(name)'." )
        end
    end
    for site in keys(sitebounds)
        if !(site in keys(siteids))
            error("Site '$(site)' provided to sitebounds does not exist on Surface '$(name)'." )
        end
        if !pyisinstance(sitebounds[site], asesf.SampleBounds)
            error("Bounds for site '$(site)' on Surface '$(name)' are not an asesurfacefinder.SampleBounds object.")
        end
    end
    # Mark unassigned coordinations and bounds as invalid placeholders.
    # These will be overwritten when put through SurfaceData.
    n_sites = length(sites)
    sitebounds = Dict{String, Any}(k => v for (k, v) in sitebounds)
    for i in 1:n_sites
        site = sites[i]
        if !(site in keys(sitecoords))
            sitecoords[site] = -1
        end
        if !(site in keys(sitebounds))
            sitebounds[site] = nothing
        end
    end

    # Extract elements.
    elements = pyconvert(Set{String}, atoms.symbols.species())

    return Surface(name, atoms, elements, sites, siteids, sitecoords, sitebounds, cutoff_mult)
end

"""
    Surface(name::String, frame::Dict{String, Any}, sitedict::Dict{String, Any}[, relative_sites=true, minimal_xy_cell=nothing, cutoff_mult=1.0, sitecoords::Dict{String, Int}=Dict(), sitebounds::Dict{String, <:Any}=Dict()])

Constructs a `Surface` compatible with ASE from an ExtXYZ `frame`.

Can take any `name`, this is just a user-facing identifier for
the surface.

`sitedict` should be a `Dict` of string-form absorption site names
bound to `[x, y]` position vectors. These can be absolute positions,
or they can be relative to the surface's unit cell. In the case of
the former, `relative_sites` should be `false` to automatically
convert to relative positions, which ASE requires for passing in
named absorption sites.

If `sitebounds` is passed, it should be a `Dict` of string-form 
absorption site names bound to instances of `asesurfacefinder.SampleBounds`.
Details can be found in ASESurfaceFinder's docs. If bounds for any
site are not provided, they will be set to ASESurfaceFinder's defaults
when passing this `Surface` to `SurfaceData`.

If `sitecoords` is passed, it should be a `Dict` of string-form
absorption site names bound to integer site coordinations. If
coordinations for any sites are not provided, they will similarly
be populated when the `Surface` is passed to `SurfaceData`.

The frame should be a minimal cell in xy, i.e it can have vertical
thickness, but should not have any periodicity in the xy plane. If
this is not the case, the minimal xy unit cell should be specified
in `minimal_xy_cell`, otherwise adsorbates will be placed incorrectly.
"""
function Surface(name::String, frame::Dict{String, Any}, sitedict::Dict; 
                 relative_sites=true, minimal_xy_cell=nothing, cutoff_mult=1.0, 
                 sitecoords::Dict{String, Int}=Dict(), sitebounds::Dict{String, <:Any}=Dict())
    atoms = frame_to_atoms(frame)
    elements = pyconvert(Set{String}, atoms.symbols.species())

    # Sanitise absorption site positions.
    for (site, sitepos) in sitedict
        if !(typeof(sitepos) <: Vector && eltype(sitepos) <: AbstractFloat && length(sitepos) == 2)
            error("Incorrectly defined absorption site '$(site)' in Surface '$(name)'.")
        end
    end

    # Convert absolute site positions into relative ones.
    if !(relative_sites)
        cell = frame["cell"][1:2, 1:2]
        invcell = inv(cell)
        for site in keys(sitedict)
            sitedict[site] = invcell * sitedict[site]
        end
    end

    atoms.info["adsorbate_info"] = pydict(; sites=pydict())
    if isnothing(minimal_xy_cell)
        # Use the full cell from the frame.
        atoms.info["adsorbate_info"]["cell"] = Py(frame["cell"][1:2, 1:2]).to_numpy()
    else
        # Use a minimal xy cell.
        atoms.info["adsorbate_info"]["cell"] = Py(minimal_xy_cell).to_numpy()
    end
    for site in keys(sitedict)
        atoms.info["adsorbate_info"]["sites"][site] = sitedict[site]
    end
    site_names = sort(pyconvert(Vector{String}, atoms.info["adsorbate_info"]["sites"].keys()))
    kinetica_sites = Dict{Int, String}(i => s for (i, s) in enumerate(site_names))
    kinetica_siteids = Dict{String, Int}(s => i for (i, s) in enumerate(site_names))
    for site in keys(sitecoords)
        if !(site in keys(kinetica_siteids))
            error("Site '$(site)' provided to sitecoords does not exist on Surface '$(name)'." )
        end
    end
    for site in keys(sitebounds)
        if !(site in keys(kinetica_siteids))
            error("Site '$(site)' provided to sitebounds does not exist on Surface '$(name)'." )
        end
        if !pyisinstance(sitebounds[site], asesf.SampleBounds)
            error("Bounds for site '$(site)' on Surface '$(name)' are not an asesurfacefinder.SampleBounds object.")
        end
    end
    # Mark unassigned coordinaitons and bounds as invalid placeholders.
    # These will be overwritten when put through SurfaceData.
    n_sites = length(sitedict)
    sitebounds = Dict{String, Any}(k => v for (k, v) in sitebounds)
    for i in 1:n_sites
        site = kinetica_sites[i]
        if !(site in keys(sitecoords))
            sitecoords[site] = -1
        end
        if !(site in keys(sitebounds))
            sitebounds[site] = nothing
        end
    end

    return Surface(name, atoms, elements, kinetica_sites, kinetica_siteids, sitecoords, sitebounds, cutoff_mult)
end


mutable struct SurfaceData
    surfaces::Vector{Surface}
    nameToSurf::Dict{String, Surface}
    nameToInt::Dict{String, Int}
    n::Int

    finder::Py
    sf_samples::Int
    sf_reject_coord::Bool
    sf_reject_bonded_h::Bool
    sf_kwargs
end

"""
    SurfaceData(surfaces::Vector{Surface}[, sf_samples=5000, reject_wrong_coordination=false, reject_bonded_hydrogens=true, kwargs...])

Constructs a container for `Surface`s and a linked ASESurfaceFinder instance.

Each `Surface` in `surfaces` can be indexed by its internal
name in `nameToSurf`, and this name can be indexed back to
a surface's integer ID with `nameToInt`.

In addition to constructing an instance of ASESurfaceFinder's
`SurfaceFinder`, constructing `SurfaceData` this way also trains
this instance to predict labels for adsorbates on high-symmetry
sites. This is currently limited to serial descriptor generation
and training only until parallelisation is stabilised. `sf_samples`
is passed through when training to specify the number of samples
that should be taken per site on each surface.

Any remaining keyword arguments are passed through to the
`SurfaceFinder` class at construction, allowing for finer control
over additional properties.

Any surface site sampling bounds or coordinations that were not
provided at `Surface` construction will be backfilled from 
the created `SurfaceFinder` instance's defaults.
"""
function SurfaceData(surfaces::Vector{Surface}; sf_samples=5000, reject_wrong_coordination=false, reject_bonded_hydrogens=true, kwargs...)
    nameToSurf = Dict(surf.name => surf for surf in surfaces)
    nameToInt = Dict(surf.name => i for (i, surf) in enumerate(surfaces))
    n = length(surfaces)

    labels = [surf.name for surf in surfaces]
    surf_atoms = [surf.atoms for surf in surfaces]
    site_coordinations = [Dict(site => c for (site, c) in surf.sitecoords if c != -1) for surf in surfaces]
    sample_bounds = [Dict(site => b for (site, b) in surf.sitebounds if !isnothing(b)) for surf in surfaces]
    full_kwargs = merge(
        Dict(
            :labels => labels, 
            :sample_bounds => sample_bounds, 
            :site_coordinations => site_coordinations,
            :verbose => false
        ),
        Dict(a => b for (a, b) in kwargs)
    )
    finder = asesf.SurfaceFinder(surf_atoms; full_kwargs...)

    # Read completed coordinations and sample bounds from SurfaceFinder.
    for i in 1:n
        n_sites = length(surfaces[i].sites)
        for j in 1:n_sites
            site = surfaces[i].sites[j]
            if surfaces[i].sitecoords[site] == -1
                surfaces[i].sitecoords[site] = pyconvert(Int, finder.site_coordinations[i-1][site])
            end
            if isnothing(surfaces[i].sitebounds[site])
                surfaces[i].sitebounds[site] = finder.sample_bounds[i-1][site]
            end
        end
    end

    # Train ASESurfaceFinder model.
    finder.train(
        samples_per_site=sf_samples,
        n_jobs=1
    )

    return SurfaceData(
        surfaces,
        nameToSurf,
        nameToInt,
        n,
        finder,
        sf_samples,
        reject_wrong_coordination,
        reject_bonded_hydrogens,
        kwargs
    )
end


"""
    add_surface!(surfdata::SurfaceData, surface::Surface)

Adds a `Surface` to `SurfaceData`.

Also adds this surface to the underlying ASESurfaceFinder
instance, invoking the retraining of its random forest
classifier.
"""
function add_surface!(surfdata::SurfaceData, surface::Surface)
    push!(surfdata.surfaces, surface)
    surfdata.n += 1
    surfdata.nameToSurf[surface.name] = surface
    surfdata.nameToInt[surface.name] = surfdata.n

    labels = [surf.name for surf in surfdata.surfaces]
    surf_atoms = [surf.atoms for surf in surfdata.surfaces]
    site_coordinations = [Dict(site => c for (site, c) in surf.sitecoords if c != -1) for surf in surfaces]
    sample_bounds = [Dict(site => b for (site, b) in surf.sitebounds if !isnothing(b)) for surf in surfaces]
    full_kwargs = merge(
        Dict(
            :labels => labels, 
            :sample_bounds => sample_bounds, 
            :site_coordinations => site_coordinations,
            :verbose => false
        ),
        Dict(a => b for (a, b) in surfdata.kwargs)
    )
    finder = asesf.SurfaceFinder(surf_atoms; full_kwargs...)

    n_sites = length(surface.sites)
    surfidx = surfdata.n
    for i in 1:n_sites
        site = surface.sites[i]
        if surface.sitecoords[site] == -1
            surfaces[surfidx].sitecoords[site] = finder.site_coordinations[surfidx-1][site]
        end
        if isnothing(surface.sitebounds[site])
            surfaces[surfidx].sitebounds[site] = finder.sample_bounds[surfidx-1][site]
        end
    end

    finder.train(
        samples_per_site=surfdata.sf_samples,
        n_jobs=1
    )
    surfdata.finder = finder
    return
end


"""
    get_surfid(smi::String)

Returns the surface ID of adsorption sites in a surface SMILES.

If a gas-phase species (i.e. a regular SMILES without a surface
tag) is passed, returns nothing.
"""
get_surfid(smi::String) = get_surfid(SpeciesStyle(smi), smi)
function get_surfid(::SurfaceSpecies, smi::String)
    m = match(r"(X\d_\d)", smi)
    return parse(Int, split(m.match[2:end], '_')[1])
end
get_surfid(::GasSpecies, smi::String) = nothing


"""
    get_surf_siteids(smi::String)

Returns the surface site IDs of adsorption sites in a surface SMILES.

If a gas-phase species (i.e. a regular SMILES without a surface
tag) is passed, returns nothing.
"""
get_surf_siteids(smi::String) = get_surf_siteids(SpeciesStyle(smi), smi)
function get_surf_siteids(::SurfaceSpecies, smi::String)
    m = 1
    offset = 0
    siteids = Int[]
    while !isnothing(m)
        m = match(r"(X\d_\d)", smi, offset+1)
        if isnothing(m) continue end
        push!(siteids, parse(Int, split(m.match[2:end], '_')[2]))
        offset = m.offset
    end
    
    return siteids
end
get_surf_siteids(::GasSpecies, smi::String) = nothing


"""
    get_surf_site_atomids(amsmi::String)
"""
get_surf_site_atomids(amsmi::String) = get_surf_site_atomids(SpeciesStyle(amsmi), amsmi)
function get_surf_site_atomids(::SurfaceSpecies, amsmi::String)
    surfid = get_surfid(amsmi)
    siteids = get_surf_siteids(amsmi)
    pt = rdChem.GetPeriodicTable()

    # Substitute out bonding to surface sites.
    site_replacements = []
    re = r"(?<=[#=])(\[X\d_\d\])" # Search for #/= before site tags
    m = match(re, amsmi)
    while !isnothing(m)
        push!(site_replacements, Pair(amsmi[m.offset-1]*m.match, m.match))
        m = match(re, amsmi, m.offset+length(m.match))
    end
    re = r"(\[X\d_\d\])(?=[#=])" # Search for #/= after site tags
    m = match(re, amsmi)
    while !isnothing(m)
        push!(site_replacements, Pair(m.match*amsmi[m.offset+length(m.match)], m.match))
        m = match(re, amsmi, m.offset+length(m.match)+1)
    end
    amsmi_subbed = replace(amsmi, site_replacements...)

    # Substitute surface site tags with unique elements.
    site_atomic_number = 100
    elem_replacements = []
    elems_used = []
    for siteid in siteids
        elem = pyconvert(String, pt.GetElementSymbol(site_atomic_number))
        push!(elem_replacements, Pair("X$(surfid)_$(siteid)", elem))
        push!(elems_used, elem)
        site_atomic_number += 1
    end
    amsmi_replaced = replace(amsmi_subbed, elem_replacements...)

    # Handle isolated hydrogen atoms manually since they cause RDKit problems.
    if amsmi_replaced == "[Fm][H:1]"
        return Dict(first(elem_replacements[1]) => 1)
    end

    mol = rdChem.MolFromSmiles(amsmi_replaced)
    elem_bonded_atomids = Dict{String, Int}()
    for atom in mol.GetAtoms()
        elem = pyconvert(String, atom.GetSymbol())
        if elem in elems_used
            neighbour = atom.GetNeighbors()[0]
            elem_bonded_atomids[elem] = pyconvert(Int, neighbour.GetAtomMapNum())
        end
    end

    elem_replacement_flipped_dict = Dict(e[2] => e[1] for e in elem_replacements)
    surftag_atomids = Dict(elem_replacement_flipped_dict[elem] => elem_bonded_atomids[elem] for elem in keys(elem_bonded_atomids))

    return surftag_atomids
end
get_surf_site_atomids(::GasSpecies, amsmi::String) = nothing

"""
    remove_surface_atoms!(frame::Dict{String, Any}, surfdata::SurfaceData, surfid::Int[, is_adsorbed::Bool=false])

Removes the atoms corresponding to the `Surface` in `surfdata.surfaces[surfid]` from `frame`.

Also removes the `cell` and `pbc` keys from `frame`. Adds
and adsorbate tag to its `info` dict if `is_adsorbed==true`,
making the resulting geometry an `AdsorbedXYZ` if this is the
case or a `FreeXYZ` if not, rather than an `OnSurfaceXYZ`.
"""
function remove_surface_atoms!(frame::Dict{String, Any}, surfdata::SurfaceData, surfid::Int, is_adsorbed::Bool=false)
    surface = surfdata.surfaces[surfid]
    elems = surface.elements
    remove_idxs = [i for (i, e) in enumerate(frame["arrays"]["species"]) if e in elems]
    keep_idxs = [i for i in 1:frame["N_atoms"] if !(i in remove_idxs)]
    for arrkey in keys(frame["arrays"])
        if frame["arrays"][arrkey] isa Vector
            deleteat!(frame["arrays"][arrkey], remove_idxs)
        else
            frame["arrays"][arrkey] = frame["arrays"][arrkey][:, keep_idxs]
        end
    end
    frame["N_atoms"] -= length(remove_idxs)
    delete!(frame, "cell")
    delete!(frame, "pbc")
    if is_adsorbed
        frame["info"]["adsorbate"] = "true"
    end
    return
end


"""
    remove_adsorbate_atoms!(frame::Dict{String, Any}, surfdata::SurfaceData, surfid::Int)

Removes the atoms corresponding to the adsorbate from `frame`.

Identifies which atoms are part of the surface from `surfdata.surfaces[surfid]`,
then removes everything else.
"""
function remove_adsorbate_atoms!(frame::Dict{String, Any}, surfdata::SurfaceData, surfid::Int)
    surface = surfdata.surfaces[surfid]
    elems = surface.elements
    keep_idxs = [i for (i, e) in enumerate(frame["arrays"]["species"]) if e in elems]
    remove_idxs = [i for i in 1:frame["N_atoms"] if !(i in keep_idxs)]
    for arrkey in keys(frame["arrays"])
        if frame["arrays"][arrkey] isa Vector
            deleteat!(frame["arrays"][arrkey], remove_idxs)
        else
            frame["arrays"][arrkey] = frame["arrays"][arrkey][:, keep_idxs]
        end
    end
    frame["N_atoms"] -= length(remove_idxs)
    return
end