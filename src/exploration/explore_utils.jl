"""
    sd, rd = import_mechanism(rdir, rcount[, max_molecularity])

Create a CRN's initial `SpeciesData` and `RxData` from a CDE generated mechanism(s).
"""
function import_mechanism(rdir::String, rcount; max_molecularity=2)
    rsmis, rxyzs, psmis, pxyzs, dHs = ingest_cde_run(rdir, rcount)
    all_smis = vcat(reduce(vcat, rsmis), reduce(vcat, psmis))
    all_xyzs = vcat(reduce(vcat, rxyzs), reduce(vcat, pxyzs))
    sd = SpeciesData(all_smis, all_xyzs)
    rd = RxData(sd, rsmis, psmis, dHs; max_molecularity=max_molecularity)
    return sd, rd
end

"""
    import_mechanism!(sd, rd, rdir, rcount[, max_molecularity])

Extend a CRN's `SpeciesData` and `RxData` from a CDE generated mechanism(s).
"""
function import_mechanism!(sd::SpeciesData, rd::RxData, rdir::String, rcount;
        max_molecularity=2)
    rsmis, rxyzs, psmis, pxyzs, dHs = ingest_cde_run(rdir, rcount)
    all_smis = vcat(reduce(vcat, rsmis), reduce(vcat, psmis))
    all_xyzs = vcat(reduce(vcat, rxyzs), reduce(vcat, pxyzs))
    push_unique!(sd, all_smis, all_xyzs)
    push!(rd, sd, rsmis, psmis, dHs; max_molecularity=max_molecularity)
    return
end