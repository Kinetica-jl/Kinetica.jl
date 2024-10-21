# Constructor for the ASE IdealGasThermo object and entropy/enthalpy
# change calculator.
# Not intended for production use, only for ensuring that Kinetica's
# dS and dH calculations match and make sense.

"""
"""
function build_ase_idealgasthermo(vib_energies::Vector{Float64}, 
                                  geometry::Int,
                                  frame::Dict{String, Any},
                                  symmetry::Int,
                                  mult::Int)
    atoms = frame_to_atoms(frame)
    geom_str = Dict(0 => "monatomic", 1 => "linear", 2 => "nonlinear")[geometry]
    spin = (mult - 1)/2

    # Atoms initialised separately to avoid cutting down number of vibrational energies.
    thermo = asethermo.IdealGasThermo(vib_energies,
                                      geom_str,
                                      potentialenergy=frame["info"]["energy_ASE"],
                                      symmetrynumber=symmetry,
                                      spin=spin)
    thermo.atoms = atoms
    return thermo
    
end

function build_ase_idealgasthermo(sd::SpeciesData, sid)
    return build_ase_idealgasthermo(sd.cache[:vib_energies][sid],
                                    sd.cache[:geometry][sid],
                                    sd.xyz[sid],
                                    sd.cache[:symmetry][sid],
                                    sd.cache[:mult][sid])
end

function build_ase_idealgasthermo(ts_cache::Dict{Symbol, Any}, rid)
    return build_ase_idealgasthermo(ts_cache[:vib_energies][rid],
                                    ts_cache[:geometry][rid],
                                    ts_cache[:xyz][rid],
                                    ts_cache[:symmetry][rid],
                                    ts_cache[:mult][rid])
end

"""
"""
function calculate_ase_idealgasthermo(calc::ASENEBCalculator, T, P)
    dS = zeros(calc.rd.nr)
    dH = zeros(calc.rd.nr)
    for rid in 1:calc.rd.nr
        S_reacs = 0.0
        H_reacs = 0.0
        mass_ts = 0.0

        for (i, sid) in enumerate(calc.rd.id_reacs[rid])
            spec_stoic = calc.rd.stoic_reacs[rid][i]
            mass_ts += spec_stoic * calc.sd.cache[:weights][sid]
            idealgasthermo = build_ase_idealgasthermo(calc.sd, sid)
            entropy = pyconvert(Float64, idealgasthermo.get_entropy(T, P, verbose=false))
            S_reacs += spec_stoic * entropy
            enthalpy = pyconvert(Float64, idealgasthermo.get_enthalpy(T, verbose=false))
            H_reacs += spec_stoic * enthalpy
        end

        idealgasthermo = build_ase_idealgasthermo(calc.ts_cache, rid)
        S_ts = pyconvert(Float64, idealgasthermo.get_entropy(T, P, verbose=false))
        H_ts = pyconvert(Float64, idealgasthermo.get_enthalpy(T, verbose=false))

        dS[rid] = S_ts - S_reacs
        dH[rid] = H_ts - H_reacs
    end

    # Convert from eV/K and eV to J/mol/K and J/mol
    dS ./= (Constants.J/Constants.mol)
    dH ./= (Constants.J/Constants.mol)

    return dS, dH
end