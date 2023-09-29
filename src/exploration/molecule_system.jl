"""
    f = spring_force(x1, x2, dx1, dx2, ks, kr, r)

Calculate spring force between two particles.
"""
function spring_force(x1, x2, dx1, dx2, ks, kd, r)
    n = norm(x1 .- x2)
    don = (x1 .- x2) / n
    return -1.0 .*(ks*(n-r) .+ kd*(dx1 .- dx2).*don) .* don
end


"""
    ddu = spring3d(du, u, k, t)

2nd order ODE around which the spring-particle solver is based.
"""
function spring3d(du, u, k, t)
    np = size(u, 1)
    ddu = zeros(Float64, np, 3)
    ks = k[end-2]
    kd = k[end-1]
    r = k[end]
    for i in 1:np
        for j in i+1:np
            u1 = u[i, :]
            u2 = u[j, :]
            du1 = du[i, :]
            du2 = du[j, :]
            force = spring_force(u1, u2, du1, du2, ks, kd, r)
            ddu[i, :] += force / k[i]
            ddu[j, :] -= force / k[j]
        end
    end
    return ddu
end


"""
    com = get_COM(mol)

Calculates center of mass (COM) of a given ExtXYZ `mol`.
"""
function get_COM(mol::Dict{String}{Any})
    com = zeros(Float64, 3)
    mass = 0.0
    for i in 1:mol["N_atoms"]
        m = atom_mass_dict[mol["arrays"]["species"][i]]
        mass += m
        com .+= (mol["arrays"]["pos"][:, i] * m)
    end
    com /= mass
    return com
end


"""
    mass = get_mass(mol)

Calculates mass of a given ExtXYZ `mol`.
"""
function get_mass(mol::Dict{String}{Any})
    mass = 0.0
    for i in 1:mol["N_atoms"]
        mass += atom_mass_dict[mol["arrays"]["species"][i]]
    end
    return mass
end


"""
    center_mols!(mols)

Centers the geometric centres of an array of molecules on the origin.
"""
function center_mols!(mols::Vector{Dict{String, Any}})
    for mol in mols
        geometric_centre = [mean(mol["arrays"]["pos"][i, :]) for i in 1:3]
        for i in 1:mol["N_atoms"]
            mol["arrays"]["pos"][:, i] -= geometric_centre
        end
    end
end


"""
    transform_mol!(mol, vec)

Transforms coordinates of `mol` by applying cartesian vector `vec`.
"""
function transform_mol!(mol::Dict{String}{Any}, vec::Vector{Float64})
    for i in 1:mol["N_atoms"]
        mol["arrays"]["pos"][:, i] += vec
    end
end


"""
    out_of_bounds = position_check(tmols)

Checks all proposed atom positions to ensure they can be read in by xTB.

xTB reads in XYZ files with a format string on positions that only allows
3 digits (including a hyphen for negative floats) before the decimal point.
This means that positions can only go from -99.999... to 100.000...
"""
function position_check(tmols::Vector{Dict{String}{Any}})
    n_mols = length(tmols)
    out_of_bounds = false
    for i in 1:n_mols
        coords = tmols[i]["arrays"]["pos"]
        if any(coords .>= 1000.0) || any(coords .<=-100.0)
            out_of_bounds = true
            break
        end
    end
    return out_of_bounds
end


"""
    too_close = proximity_check(tmols, dmin)

Checks all atoms in each molecule in `tmols` are at least 
`dmin` Angstroms apart from atoms in other molecules.
"""
function proximity_check(tmols::Vector{Dict{String}{Any}}, dmin::Float64)
    n_mols = length(tmols)
    too_close = false
    for i in 1:n_mols, j in i+1:n_mols
        ni = tmols[i]["N_atoms"]
        nj = tmols[j]["N_atoms"]
        prox = zeros(Float64, ni, nj)
        for k in 1:ni, l in 1:nj
            prox[k, l] = norm(tmols[i]["arrays"]["pos"][:, k] - tmols[j]["arrays"]["pos"][:, l])
        end
        if !all(prox .>= dmin)
            too_close = true
            break
        end
    end
    return too_close
end


"""
    tmols = molsys_opt(mols, dmin, maxiters)

Optimises positions of molecules in `mols` to ensure they are all at least `dmin` Angstroms apart.

Creates and solves an N-body spring-driven particle system and
transforms molecular coordinates to these particles to check
for proximity.
"""
function molsys_opt(mols::Vector{Dict{String}{Any}}, dmin::Float64, maxiters::Int)
    @debug "Starting new molsys optimisation"
    np = length(mols)

    x0 = zeros(Float64, np, 3) .+ reshape(rand(np*3), np, 3)
    dx0 = zeros(Float64, np, 3)
    k = [get_mass(mol) for mol in mols]
    k = vcat(k, [
        2.0,    # Spring constant
        0.75,   # Damping constant
        40.0    # Rest length
    ])
    tspan = (0.0, 1e5)

    condition(u, t, integrator) = norm(u[1:np*3]) <= 1e-4 ? true : false
    affect!(integrator) = terminate!(integrator)
    callback = DiscreteCallback(condition, affect!)
    prob = SecondOrderODEProblem(spring3d, dx0, x0, tspan, k, callback=callback)

    too_close = true
    oob = true
    counter = 0
    r_adj = 0.0
    tmols = []
    while too_close || oob
        counter += 1
        if counter >= maxiters
            error("Max iterations exceeded in molsys_opt().")
        end

        new_k = copy(k)
        new_k[end] += r_adj
        new_u0 = copy(prob.u0)
        new_u0[(np*3)+1:end] = zeros(Float64, np, 3) .+ reshape(rand(np*3), np, 3)
        new_prob = remake(prob; u0=new_u0, p=new_k)
        @debug "Iter $counter - Spring rest length = $(new_prob.p[end])"
        sol = solve(new_prob, Tsit5())

        solmat = reduce(vcat, sol.u')
        xmat = solmat[:, (np*3)+1:end]
        px = [xmat[end, 3(i-1)+1:3(i-1)+3] for i in 1:np]

        tmols = deepcopy(mols)
        for i in 1:np
            transform_mol!(tmols[i], px[i])
        end
        oob = position_check(tmols)
        too_close = proximity_check(tmols, dmin)
        if oob && too_close
            @debug "Transformed molecules too close and out of bounds"
            r_adj -= 5.0
        elseif oob && !too_close
            @debug "Transformed molecules out of bounds"
            r_adj -= 10.0
        elseif too_close && !oob
            @debug "Transformed molecules too close"
            r_adj += 10.0
        end
    end
    @debug "Finished molsys optimisation"
    return tmols
end


"""
    mol = combine_mols(tmols)

Combines a vector of ExtXYZ Atoms dicts into a single Atoms Dict.

Copies each constituent XYS into the new unified XYZ blindly, so
resoponsibility is on the user to check there is no overlap in any
of the coordinates.
"""
function combine_mols(tmols::Vector{Dict{String}{Any}})
    nmols = length(tmols)
    newmol = copy(tmols[1])
    for i in 2:nmols
        newmol["N_atoms"] += tmols[i]["N_atoms"]
        newmol["arrays"]["pos"] = hcat(newmol["arrays"]["pos"], tmols[i]["arrays"]["pos"])
        newmol["arrays"]["species"] = vcat(newmol["arrays"]["species"], tmols[i]["arrays"]["species"])
    end
    return newmol
end


"""
    system_from_smiles(smiles[, saveto, dmin])

Forms a single XYZ system out of the molecules in `smiles`.

Useful for making unified molecular systems with no overlap
for feeding into CDE. `dmin` represents the minimum 
molecule-molecule distance that should be allowed.

If the argument `saveto` is provided, outputs the optimised
system to a file at this path. If not, returns the optimised
system as a single ExtXYZ dict.
"""
function system_from_smiles(smiles::Vector{String}, 
        saveto::String; dmin::Float64=5.0, maxiters::Int=200)

    mol = system_from_smiles(smiles; dmin=dmin, maxiters=maxiters)
    write_frame(saveto, mol)
end

function system_from_smiles(smiles::Vector{String}; 
        dmin::Float64=5.0, maxiters::Int=200)

    for (i, smi) in enumerate(smiles)
        xyz_from_smiles(smi, joinpath(dirname(saveto), "tmp_$i.xyz"))
    end
    mols = [read_frame(joinpath(dirname(saveto), "tmp_$i.xyz")) for i in 1:length(smiles)]
    for i in 1:length(smiles)
        rm(joinpath(dirname(saveto), "tmp_$i.xyz"))
    end

    center_mols!(mols)
    tmols = molsys_opt(mols, dmin, maxiters)
    mol = combine_mols(tmols)
    return mol
end


"""
"""
function system_from_mols(mols::Vector{Dict{String}{Any}}, saveto::String; 
        dmin::Float64=5.0, maxiters::Int=200)

    mol = system_from_mols(mols; dmin=dmin, maxiters=maxiters)
    write_frame(saveto, mol)
end

function system_from_mols(mols::Vector{Dict{String}{Any}}; 
        dmin::Float64=5.0, maxiters::Int=200)

    center_mols!(mols)
    tmols = molsys_opt(mols, dmin, maxiters)
    mol = combine_mols(tmols)
    return mol
end

atom_mass_dict = Dict{String, Float64}(
    "H" => 1.008,
    "He" => 4.002602,
    "Li" => 6.94,
    "Be" => 9.0121831,
    "B" => 10.81,
    "C" => 12.011,
    "N" => 14.007,
    "O" => 15.999,
    "F" => 18.99840316,
    "Ne" => 20.1797
)