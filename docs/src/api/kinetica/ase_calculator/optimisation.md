# Kinetica.jl API - ASE Interface

## Optimisation and Property Calculation

### Basic Species Properties

```@docs
Kinetica.get_mult
Kinetica.get_mult!
Kinetica.get_charge
Kinetica.get_charge!
Kinetica.get_formal_charges
Kinetica.get_formal_charges!
Kinetica.get_initial_magmoms
Kinetica.get_initial_magmoms!
Kinetica.correct_magmoms_for_mult!
Kinetica.get_hydrogen_idxs
```

### Geometry Optimisation

```@docs
Kinetica.geomopt!
Kinetica.kabsch_fit!
Kinetica.permute_hydrogens!
```

### NEB

```@docs
Kinetica.get_initial_sys_mult
Kinetica.get_rxn_mult
Kinetica.neb
Kinetica.highest_energy_frame
```

### Vibrational Analysis

```@docs
Kinetica.calc_species_vibrations!
Kinetica.calc_ts_vibrations!
```