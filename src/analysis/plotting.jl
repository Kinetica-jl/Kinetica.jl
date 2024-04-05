@recipe function f(res::ODESolveOutput;
        tunit="s",
        label_above=0.1,
        ignore_species=nothing,
        ignore_below=nothing
    )

    xlabel --> "Time / $tunit"
    ylabel --> "Concentration / mol dm⁻³"
    legend_position --> :topright
    palette --> :default

    umat = reduce(vcat, res.sol.u')
    x = res.sol.t
    colour_count = 1
    for i in axes(umat, 2)
        if typeof(ignore_species) <: Vector && res.sd.toStr[i] in ignore_species
            continue
        end

        y = umat[:, i]
        if !isnothing(ignore_below) && maximum(y) < ignore_below
            continue
        end
        
        important_series = maximum(y) >= label_above
        slabel = important_series ? res.sd.toStr[i] : ""
        scolor = important_series ? colour_count : "grey"
        if important_series colour_count += 1 end
        coords = [(ix, iy) for (ix, iy) in zip(x, y)]
        
        @series begin
            seriestype := :path
            label --> slabel
            color --> scolor
            coords
        end
    end

    primary := false
    ()
end


@userplot ConditionsPlot
@recipe function f(cp::ConditionsPlot; tunit="s")
    res = cp.args[1]
    if typeof(res) <: ConditionSet
        cs = res
    elseif typeof(res) <: ODESolveOutput
        cs = res.conditions
    else
        error("First argument of conditionsplot must be a ConditionSet or an ODESolveOutput containing a ConditionSet.")
    end

    csym = cp.args[2]
    profile = get_profile(cs, csym)
    if !isvariable(profile)
        error("Profile for condition $csym is not variable.")
    end

    xlabel --> "Time / $tunit"
    ylabel --> csym
    legend --> false
    x := profile.sol.t
    y := reduce(vcat, profile.sol.u)

    ()
end


function sort_species_final(res::ODESolveOutput)
    u_selected = res.sol.u[end]
    u_sorted = sort(u_selected, rev=true)
    ids_sorted = sortperm(u_selected, rev=true)

    return ids_sorted, u_sorted
end

@userplot FinalConcPlot
@recipe function f(fcp::FinalConcPlot; quantity=:conc, n_top=10, highlight_radicals=false, ignore_species=nothing, xscale=:identity)
    res = fcp.args[1]
    if !(typeof(res) <: ODESolveOutput)
        error("First argument of finalconcplot must be an instance of ODESolveOutput.")
    end
    if !(quantity in [:percent, :conc])
        error("Unknown value for `quantity`, must be one of [:percent, :conc].")
    end

    sorted_ids, sorted_concs = sort_species_final(res)
    if quantity == :percent
        sorted_quantity = sorted_concs / sum(sorted_concs) * 100
    else
        sorted_quantity = sorted_concs
    end
    sorted_names = [res.sd.toStr[i] for i in sorted_ids]

    top_names = sorted_names[1:n_top]
    top_quantity = sorted_quantity[1:n_top]
    others_quantity = sum(sorted_quantity[n_top+1:end])

    if !isnothing(ignore_species)
        for spec in ignore_species
            if spec in top_names
                remid = findfirst(x -> x == spec, top_names)
                deleteat!(top_names, remid)
                deleteat!(top_quantity, remid)
                push!(top_names, sorted_names[n_top+1])
                push!(top_quantity, sorted_quantity[n_top+1])
                others_quantity -= sorted_quantity[n_top+1]
            end
        end
    end

    push!(top_names, "Others")
    push!(top_quantity, others_quantity)
    yvals = collect(1:length(top_names))

    colours = [1 for _ in 1:length(top_names)]
    if highlight_radicals
        for i in axes(colours)[1]
            if occursin("[", top_names[i]) && top_names[i] != "[H][H]"
                bracketpos = Vector{Int}()
                valid_radical = false
                for (j, cha) in enumerate(top_names[i])
                    if cha == '['
                        push!(bracketpos, j)
                    end
                end
                for pos in bracketpos
                    if !(top_names[i][pos+1:pos+2] == "C@")
                        valid_radical = true
                    end
                end
                if valid_radical
                    colours[i] = 2
                end
            end
        end
    end

    seriestype := :bar
    orientation := :h
    xscale := xscale
    y := top_quantity
    x := yvals
    yticks := (yvals, top_names)
    yflip := true
    xlabel := "Concentration / mol dm⁻³"
    legend := false
    color := colours

    if xscale == :identity
        xlims := (0.0, top_quantity[1])
    else
        # Deal with negative numbers due to numerical noise.
        min_top_quantity = maximum(top_quantity)
        for q in top_quantity
            if 0 < q < min_top_quantity
                min_top_quantity = q
            end
        end
        low_lim = 10^floor(log(10, min_top_quantity))
        xlims := (low_lim, 10^ceil(log(10, top_quantity[1])))
        fillrange := low_lim
    end

    ()
end