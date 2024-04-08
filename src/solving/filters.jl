mutable struct RxFilter
    filters::Vector
    keep_filtered::Bool
end

"""
    RxFilter(filters[, keep_filtered=false])

Data container for CRN filters.

Defines a set of functions which can be used on a network to construct
a set of masks of reactions. These reactions can then be excluded from
a network, or their inverse can be taken to exclude all other reactions
form a network.

Contains fields for:
* Array of filter functions, each taking a tuple of `::SpeciesData` and `::RxData` as arguments (`filters`)
* Whether to remove or keep masked reactions (`keep_filtered`)

Can be constructed blank (`rf = RxFilter()`) to obtain a mask of all
reactions, which are then kept. Can also be constructed as 
`rf = RxFilter(filters)` to default to removing the filtered reactions.
"""
function RxFilter(filters::Vector; keep_filtered=false)
    return RxFilter(filters, keep_filtered)
end

function RxFilter()
    return RxFilter([(sd, rd) -> [false for _ in 1:rd.nr]], false)
end


"""
    get_filter_mask(rf::RxFilter, sd::SpeciesData, rd::RxData)

Calculates a combined reaction mask from all filters in `rf`.

Inverts the combined mask if `rf.keep_filtered == true`.
"""
function get_filter_mask(rf::RxFilter, sd::SpeciesData, rd::RxData)
    if length(rf.filters) == 0 error("RxFilter has not filter functions defined.") end
    
    inv_mask = .~rf.filters[1](sd, rd)
    for i in 2:length(rf.filters)
        filter_mask = rf.filters[i](sd, rd)
        inv_mask .*= .~filter_mask
    end
    mask = .~inv_mask

    if rf.keep_filtered mask = .~mask end
    return mask
end