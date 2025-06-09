"""
    calc_species_vibrations(sd::SpeciesData, sid, calc_builder[, calcdir="./", refresh=false, delta=0.01, ivetol=0.1, kwargs...])

Calculates vibrational energies of a species.

Given a species with an optimised geometry in `sd.xyz[sid]`,
uses ASE to calculate vibrational energies with a finite
difference approximation of its Hessian, using the ASE
calculator constructed through `calc_builder`.

Writes the vibrational energies in eV as an array into 
`sd.cache[:vib_energies]`. Defaults to not overwriting 
already calculated vibrational energies for a species,
instead returnong immediately. If recalculations are
desired, use `refresh=true`. 

If a species has imaginary frequencies, these can be
ignored below a vibrational energy threshold `ivetol`,
which removes them from the final vibrational energy
array. All imaginary frequencies can be ignored by
passing `ivetol=0.0`.
"""
function calc_species_vibrations!(sd::SpeciesData, sid, calc_builder; calcdir::String="./", refresh=false, delta=0.01, ivetol=0.1, kwargs...)
    if sid in keys(sd.cache[:vib_energies]) && !(refresh)
        @debug "Species $sid has vibrations cached, skipping."
        return
    end

    if sd.cache[:geometry][sid] == 0
        @debug "Species $sid is monoatomic, skipping vibrational analysis."
        sd.cache[:vib_energies][sid] = []
        return
    end

    atoms = frame_to_atoms(sd.xyz[sid], sd.cache[:formal_charges][sid], sd.cache[:initial_magmoms][sid])
    atoms.calc = calc_builder(calcdir, sd.cache[:mult][sid], sd.cache[:charge][sid], kwargs...)

    vibdir = joinpath(calcdir, "vib")
    if !(isdir(vibdir)) mkdir(vibdir) end
    vib = asevib.Vibrations(atoms, delta=delta)
    vib.run()
    vib_energies = vib.get_energies()
    @debug "vib_energies = $(vib_energies)"
    if SpeciesStyle(sd.toStr[sid]) isa SurfaceSpecies
        # Do nothing, we want all vibrational energies for harmonic limit TST
    elseif sd.cache[:geometry][sid] == 1
        vib_energies = vib_energies[pyslice(-(3*sd.xyz[sid]["N_atoms"] - 5), pylen(vib_energies))]
    elseif sd.cache[:geometry][sid] == 2
        vib_energies = vib_energies[pyslice(-(3*sd.xyz[sid]["N_atoms"] - 6), pylen(vib_energies))]
    else
        throw(ErrorException("Unknown geometry for species $sid: $(sd.cache[:geometry][sid])"))
    end
    jlve = pyconvert(Vector{ComplexF64}, vib_energies)
    @debug "vib_energies (reduced) = $jlve"

    if ivetol <= 0.0
        n_ve = length(jlve)
        jlve = [z.re for z in jlve if z.re > 0.0]
        n_ive = n_ve - length(jlve)
        if n_ive > 0 @debug "Removed $(n_ive) imaginary modes." end
    else
        if any([z.im > ivetol for z in jlve])
            throw(ErrorException("Imaginary frequency detected in geometry of species $sid."))
        else
            jlve = [z.re for z in jlve if z.re > 0.0]
        end
    end
    @debug "vib_energies (processed) = $jlve"

    sd.cache[:vib_energies][sid] = jlve
    rm(vibdir, recursive=true)
    return
end

"""
    calc_ts_vibrations(ts_cache::Dict{Symbol, Any}, rid, calc_builder[, calcdir="./", delta=0.01, ivetol=0.1, kwargs...])

Calculates vibrational energies of a transition state.

Given a TS with an optimised geometry in `ts_cache[:xyz][rid]`,
uses ASE to calculate vibrational energies with a finite 
difference apporoximation of its Hessian, using the ASE
calculator constructed through `calc_builder`.

Writes the vibrational energies in eV as an array into 
`ts_cache[:vib_energies]`.

If a species has imaginary frequencies, these can be
ignored below a vibrational energy threshold `ivetol`,
which removes them from the final vibrational energy
array. All imaginary frequencies can be ignored by
passing `ivetol=0.0`.
"""
function calc_ts_vibrations!(ts_cache::Dict{Symbol, Any}, rid, calc_builder; calcdir::String="./", delta=0.01, ivetol=0.1, kwargs...)
    atoms = frame_to_atoms(ts_cache[:xyz][rid], ts_cache[:xyz][rid]["info"]["formal_charges"], ts_cache[:xyz][rid]["info"]["initial_magmoms"])
    atoms.calc = calc_builder(calcdir, ts_cache[:mult][rid], ts_cache[:charge][rid], kwargs...)

    vibdir = joinpath(calcdir, "vib")
    if !(isdir(vibdir)) mkdir(vibdir) end
    vib = asevib.Vibrations(atoms, delta=delta)
    vib.run()
    vib_energies = vib.get_energies()
    @debug "vib_energies = $(vib_energies)"
    if !isempty(ts_cache[:ads_xyz][rid])
        # Do nothing, we want all vibrational energies for harmonic limit TST
    elseif ts_cache[:geometry][rid] == 1
        vib_energies = vib_energies[pyslice(-(3*ts_cache[:xyz][rid]["N_atoms"] - 5), pylen(vib_energies))]
    elseif ts_cache[:geometry][rid] == 2
        vib_energies = vib_energies[pyslice(-(3*ts_cache[:xyz][rid]["N_atoms"] - 6), pylen(vib_energies))]
    else
        throw(ErrorException("Unknown geometry for TS $rid: $(ts_cache[:geometry][rid])"))
    end
    jlve = pyconvert(Vector{ComplexF64}, vib_energies)
    @debug "vib_energies (reduced) = $jlve"
    
    if ivetol <= 0.0
        n_ve = length(jlve)
        jlve = [z.re for z in jlve if z.re > 0.0]
        n_ive = n_ve - length(jlve)
        if n_ive > 0 @debug "Removed $(n_ive) imaginary modes." end
    else
        if any([z.im > ivetol for z in jlve])
            throw(ErrorException("Imaginary frequency detected in geometry of TS $rid."))
        else
            jlve = [z.re for z in jlve if z.re > 0.0]
        end
    end
    @debug "vib_energies (processed) = $jlve"

    push!(ts_cache[:vib_energies], jlve)
    rm(vibdir, recursive=true)
    return
end