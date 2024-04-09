# Kinetica.jl API

## Conditions

### `ConditionSet`

```@docs
ConditionSet
get_profile(::ConditionSet, ::Symbol)
Kinetica.get_initial_conditions(::ConditionSet)
isstatic
isvariable
get_tstops(::ConditionSet)
get_t_final(::ConditionSet)
solve_variable_conditions!
```

### Static Condition Profiles

```@docs
Kinetica.StaticConditionProfile
```

### Variable Condition Profiles

```@docs
Kinetica.create_discrete_tstops!
```

### Directly Variable Condition Profiles

```@docs
Kinetica.solve_variable_condition!(::Kinetica.AbstractDirectProfile, ::ODESimulationParams)
NullDirectProfile
NullDirectProfile()
LinearDirectProfile
LinearDirectProfile()
```

### Gradient-Variable Condition Profiles

```@docs
Kinetica.solve_variable_condition!(::Kinetica.AbstractGradientProfile, ::ODESimulationParams)
LinearGradientProfile
LinearGradientProfile()
DoubleRampGradientProfile
DoubleRampGradientProfile()
```