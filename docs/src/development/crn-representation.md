# [CRN Representation](@id crn_representation_page)

!!! note "In Progress"
    This page is still under construction. Check back later!

When adding new reactions to a CRN, Kinetica performs a few checks. Reactions are discarded if they:

* Only contain conformational changes. Species are internally represented by SMILES, so any changes not representable by SMILES (e.g. single bond rotations) lead to species being classed as the same. Note that E/Z isomers and enantiomers are distinct from one another in SMILES.
* Exceed the maximum molecularity. By default this is 2, so only unimolecular and bimolecular reactions are allowed to enter a CRN. This is controllable by the `max_molecularity` argument when constructing or pushing to [`RxData`](@ref).
* Already exist in the [`RxData`](@ref) they are being added to. This is checked by comparing a reaction's unique hash to those currently in the CRN.

<!-- Need to mention how reverse reactions are automatically added during CDE ingest somewhere. -->