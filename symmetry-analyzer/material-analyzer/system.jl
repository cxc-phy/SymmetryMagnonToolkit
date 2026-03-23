include("../math-tool/mathtool.jl")
include("../symmetry-operation-tool/symmetry.jl")
include("../representation-tool/realspace.jl")
include("../representation-tool/kspace.jl")
include("../representation-tool/irrepresentation.jl")
include("../representation-tool/grouprepresentation.jl")

######################################################################################################################
###################################################SystemStructure####################################################
######################################################################################################################
module System

###################################################SystemStructure####################################################
######################################################################################################################
module SystemStructure
using LinearAlgebra
using ...MathTool

export get_structure
export read_INFO_for_FORTRAN

function get_structure(INFOname::String, POSCARname::String)::Dict{String, Any}
    cell = Dict{String, Any}()
    read_POSCAR!(cell, POSCARname)
    read_INFO!(cell, INFOname)
    return cell
end

function read_POSCAR!(cell, POSCARname::String)
    file = open(POSCARname, "r")
    content = readlines(file)
    lattice = zeros(Float64, 3, 3)
    atoms = Vector{String}()
    coords = Vector{Vector{Float64}}()
    scale = 0
    species = []
    for line in eachindex(content)
        if line == 2
            scale = parse(Float64, content[line])
        elseif line in 3:5
            basis = parse.(Float64, split(content[line]))
            lattice[:, line-2] = scale * basis[:]
        elseif line == 6
            append!(species, split(content[line]))
        elseif line == 7
            nums = parse.(Int64, split(content[line]))
            for i in eachindex(species)
                append!(atoms, repeat([species[i]], nums[i]))
            end
        elseif line >= 9
            push!(coords, round.(parse.(Float64, split(content[line])[1:3]), digits=5))
        end
    end
    close(file)

    if !(length(atoms) == length(atoms) == length(coords))
        error("Error occurs at read_POSCAR!")
    end
    cell["lattice"] = lattice
    cell["atom"] = atoms
    cell["frac_coord"] = coords
    cell["cartesian_coord"] = [round.(lattice*coord, digits=5) for coord in coords]
    cell["size"] = length(atoms)
end

"""
mag: 0 false; 1 lattice vector notation; 2 cartesian vector notation
"""
function read_INFO!(cell, INFOname::String)
    file = open(INFOname, "r")
    content = readlines(file)
    mag = nothing
    soc = nothing
    readstats = Dict("magmom" => false)
    magmoms = []
    for i in eachindex(content)
        #= processing line=#
        n = findfirst('#', content[i])
        if n == nothing
            line = strip(content[i])
        else
            line = strip(content[i][1:n-1])
        end
        if line == "" continue end
        #= reading info=#
        line = strip.(split(line, "="))
        if contains(line[1], "MAGNETIC")
            mag = parse(Bool, line[2])
        elseif contains(line[1], "SOC")
            soc = parse(Bool, line[2])
        elseif contains(line[1], "MAGMOM")
            if mag == 0 continue end
            for key in eachindex(readstats)
                if key == "magmom"
                    readstats[key] = true
                else
                    readstats[key] = false
                end
            end
            magmom = split(line[2])
            @goto parse_magmom
        else
            if readstats["magmom"]
                magmom = split(line[1])
                @label parse_magmom
                for i in eachindex(magmom)
                    #println(magmom[i])
                    if magmom[i] == "*"
                        multiplier = parse(Int64, magmom[i-1])
                        numstr = magmom[i+1]
                        magmom = [magmom[1:i-2]; repeat([numstr], multiplier); magmom[i+2:end]]
                    elseif contains(magmom[i], "*")
                        multiplier = parse(Int64, split(magmom[i], "*")[1])
                        numstr = split(magmom[i], "*")[2]
                        magmom = [magmom[1:i-1]; repeat([numstr], multiplier); magmom[i+1:end]]
                    end
                end

                append!(magmoms, parse.(Float64, magmom))
            end
        end
    end
    close(file)
    magmoms = [[magmoms[j] for j = 3*i-2:3*i] for i = 1:Int64(length(magmoms)/3)]
    #lattice_magmoms = [magmom/norm(magmom) for magmom in magmoms]

    cell["soc"] = soc
    if mag != 0
        cell["mag"] = true
        if length(magmoms) != cell["size"]
            error("Error occurs at read_INFO!: number of input magmoms in INFO is not compatible with cell size")
        end
        if mag == 1
            cell["magmom"] = magmoms
        elseif mag == 2
            cell["cartesian_magmoms"] = [magmom/norm(magmom) for magmom in magmoms]
        end
        #cell["cartesian_magmoms"] = [round.(cell["lattice"]*magmom, digits=5)/norm(cell["lattice"]*magmom) for magmom in lattice_magmoms]
        cell["magmom_arrangement"] = get_magmom_arrangement(magmoms)
        #rots = get_rotation_matrix_from_global_axis_to_local_axis(cell)
        #cell["rots"] = rots
    else
        cell["mag"] = false
    end
end

"""
get_magmom_arrangement(magmom) -> Int64
0:nonmagnetic
1:collinearFM
2:collinearAFM
3:coplanar
4:noncoplanar
"""
function get_magmom_arrangement(magmoms; tor=1e-3)::Int64
    magmoms = [magmom for magmom in magmoms if any(abs.(magmom) .> tor)]
    arrgmt = 0
    if magmoms == [] return arrgmt end
    arrgmt += 1
    dir_ori = [0.0, 0.0, 0.0]
    for magmom in magmoms
        if all(abs.(dir_ori) .< tor)
            dir_ori = magmom
        else
            dir = magmom
            if arrgmt == 1
                if is_vector_paralleled(dir_ori, dir) continue end
                if is_vector_antiparalleled(dir_ori, dir)
                    arrgmt += 1
                else
                    arrgmt += 2
                    dir_ori = cross(dir_ori, dir)
                end
            elseif arrgmt == 2
                if is_vector_paralleled(dir_ori, dir) || is_vector_antiparalleled(dir_ori, dir) continue end
                arrgmt += 1
                dir_ori = cross(dir_ori, dir)
            elseif arrgmt == 3
                if is_vector_perpendicular(dir_ori, dir) continue end
                arrgmt +=1
            else
                break
            end
        end
    end
    return arrgmt
end

"""
read_INFO_for_FORTRAN
"""
function read_INFO_for_FORTRAN(INFOname)
    file = open(INFOname, "r")
    content = readlines(file)
    Jex = nothing
    sym = nothing
    nbr_order = 0
    Jorder = 0
    Jineq = 0
    Jvalue = []
    for i in eachindex(content)
        #= processing line=#
        n = findfirst('#', content[i])
        if n == nothing
            line = strip(content[i])
        else
            line = strip(content[i][1:n-1])
        end
        if line == "" continue end
        #= reading info=#
        line = strip.(split(line, "="))
        if contains(line[1], "J_value")
            Jex = true
            if contains(line[2], "*")
                Jorder, Jineq = parse.(Int64, split(line[2], "*"))
            else
                Jorder = parse(Int64, line[2])
            end
        elseif contains(line[1], "NBR_ORDER")
            if contains(line[2], "*") line[2] = split(line[2], *)[1] end
            nbr_order = parse(Int64, line[2])
        elseif contains(line[1], "SYM")
            sym = parse.(Bool, replace(line[2], "."=>""))
        else
            if Jex == true && length(Jvalue) != Jorder
                push!(Jvalue, parse.(Float64, split(line[1])))
            end
        end
    end
    infofort = Dict("Jorder"=>Jorder, "Jvalue"=>Jvalue, "nbr_order"=>nbr_order, "sym"=>sym)
    return infofort
end

"""
get_rotation_matrix_from_global_axis_to_local_axis(cell; tor=1e-5) -> Vector{Matrix{Float64}}
"""
function get_rotation_matrix_from_local_axis_to_global_axis(cell; tor=1e-5)::Vector{Matrix{Float64}}
    rots = Vector{Matrix{Float64}}()
    for magmom in cell["cartesian_magmoms"]
        r = sqrt(sum(magmom.^2))
        d = sqrt(sum(magmom[1:2].^2))
        theta = acos(magmom[3]/r)
        if d < tor
            phi = 0.0
        elseif magmom[2] > -tor
            phi = acos(magmom[1]/d)
        else
            phi = 2*pi-acos(magmom[1]/d)
        end
        rot = [
             cos(phi)*cos(theta) -sin(phi) cos(phi)*sin(theta)
             sin(phi)*cos(theta)  cos(phi) sin(phi)*sin(theta)
            -sin(theta)           0.0      cos(theta)
        ]
        push!(rots, round.(rot, digits=5))
    end
    return rots
end

end# (module SystemStructure)

#####################################################SystemSite#######################################################
######################################################################################################################
module SystemSite
using LinearAlgebra
using ...MathTool
using ...Symmetry: spinspasym, spinptsym, spasym, ptsym

export get_cartesian_vector_between_sites
export is_site_equivalent
export get_environment_of_site
export is_environment_same
export perform_operation_on_site
export get_full_sites
export get_full_magnetic_sites
export get_site_on_nth_cell
export get_index_and_cell_of_a_site
export get_a_site
export is_site_same
export is_atom_same
export is_coord_same
export is_magmom_same

#=
get_cartesian_vector_between_sites(cell, i::Int64, j::Int64, n::Vector{Int64}) -> Vector{Float64}
is_site_equivalent(i, j, cell) -> Bool
get_environment_of_site(i, cell) -> Dict{String, Any}
is_environment_same(env1, env2; mag::Bool) -> Bool
perform_operation_on_site(cell, sym, (i, n)::Tuple{Int64, Vector{Int64}}) -> Tuple{Int64, Vector{Int64}}
perform_operation_on_site(sym::spinspasym, site)
perform_operation_on_site(sym::spinptsym, site)
perform_operation_on_site(sym::spasym, site)
perform_operation_on_site(sym::ptsym, site)
get_full_sites(cell; mag::Bool)
get_site_on_nth_cell(cell, i::Int64, n::Vector{Int64}; cartesian::Bool=false, mag::Bool=cell["mag"])
get_index_and_cell_of_a_site(cell, site; cartesian::Bool=false) -> Tuple{Int64, Vector{Int64}}
get_a_site(i, cell; mag::Bool)
is_site_same(site1, site2) -> Bool
is_atom_same(atom1, atom2) -> Bool
is_coord_same(coord1, coord2; tor=1e-3) -> Bool
is_magmom_same(magmom1, magmom2; tor=1e-3) -> Bool
=#


"""
get_cartesian_vector_between_sites(cell, i::Int64, j::Int64, n::Vector{Int64}) -> Vector{Float64}
"""
function get_cartesian_vector_between_sites(cell, i::Int64, j::Int64, n::Vector{Int64})::Vector{Float64}
    sitei = get_site_on_nth_cell(cell, i, [0, 0, 0], cartesian=true)
    sitej = get_site_on_nth_cell(cell, j, n, cartesian=true)
    cartesian_vector = sitej[2]-sitei[2]
    return cartesian_vector
end

"""
is_site_equivalent(i, j, cell) -> Bool
"""
function is_site_equivalent(i, j, cell)::Bool
    if cell["mag"]
        atomi = cell["atom"][i]
        atomj = cell["atom"][j]
        magmomi = cel["magmom"][i]
        magmomj = cel["magmom"][j]
        if !(is_atom_same(atomi, atomj) && is_magmom_same(magmomi, magmomj)) return false end
    else
        atomi = cell["atom"][i]
        atomj = cell["atom"][j]
        if !is_atom_same(atomi, atomj) return false end
    end
    envi = get_environment_of_site(i, cell)
    envj = get_environment_of_site(j, cell)
    if is_environment_same(envi, envj, mag=false)
        return true
    else
        return false
    end
end

"""
get_environment_of_site(i, cell) -> Dict{String, Any}
"""
function get_environment_of_site(i, cell)::Dict{String, Any}
    env = Dict{String, Any}()
    if cell["mag"]
        for j in 1:cell["size"]
            atom = cell["atom"][j]
            if !haskey(env, atom) env[atom]=[] end
            vec = mod.(round.(cell["frac_coord"][j]-cell["frac_coord"][i], digits=5), 1)
            magmom = cell["magmom"][j]
            push!(env[atom], (r=vec, m=magmom))
        end
    else
        for j in 1:cell["size"]
            atom = cell["atom"][j]
            if !haskey(env, atom) env[atom]=[] end
            vec = mod.(round.(cell["frac_coord"][j]-cell["frac_coord"][i], digits=5), 1)
            push!(env[atom], vec)
        end
    end
    return env
end

"""
is_environment_same(env1, env2; mag::Bool) -> Bool
"""
function is_environment_same(env1, env2; mag::Bool)::Bool
    if mag
        atoms = collect(keys(env1))
        for atom in atoms
            if !haskey(env2, atom) return false end
            for (coord1, magmom1) in env1[atom]
                for (coord2, magmom2) in env2[atom]
                    if is_coord_same(coord1, coord2) && is_magmom_same(magmom1, magmom2) @goto next1 end
                end
                return false
                @label next1
            end
        end
    else
        atoms = collect(keys(env1))
        for atom in atoms
            if !haskey(env2, atom) return false end
            for coord1 in env1[atom]
                for coord2 in env2[atom]
                    if is_coord_same(coord1, coord2) @goto next2 end
                end
                return false
                @label next2
            end
        end
    end
    return true
end

"""
perform_operation_on_site(cell, sym, (i, n)::Tuple{Int64, Vector{Int64}}) -> Tuple{Int64, Vector{Int64}}
perform_operation_on_site(sym::spinspasym, site)
perform_operation_on_site(sym::spinptsym, site)
perform_operation_on_site(sym::spasym, site)
perform_operation_on_site(sym::ptsym, site)
"""
function perform_operation_on_site(cell, sym, (i, n)::Tuple{Int64, Vector{Int64}})::Tuple{Int64, Vector{Int64}}
    site = get_site_on_nth_cell(cell, i, n)
    site_rot = perform_operation_on_site(sym, site)
    return get_index_and_cell_of_a_site(cell, site_rot)
end
function perform_operation_on_site(sym::spinspasym, site)
    site_rot = (site[1], round.(sym.r*site[2]+sym.t, digits=5), sym.s*site[3])
    return site_rot
end
function perform_operation_on_site(sym::spinptsym, site)
    site_rot = (site[1], round.(sym.r*site[2], digits=5), sym.s*site[3])
    return site_rot
end
function perform_operation_on_site(sym::spasym, site)
    site_rot = (site[1], round.(sym.r*site[2]+sym.t, digits=5))
    return site_rot
end
function perform_operation_on_site(sym::ptsym, site)
    site_rot = (site[1], sym*site[2])
    return site_rot
end

"""
get_full_sites(cell; cartesian::Bool=false, mag::Bool=cell["mag"])
"""
function get_full_sites(cell; cartesian::Bool=false, mag::Bool=cell["mag"])
    if mag
        if cartesian
            sites = [(cell["atom"][i], cell["cartesian_coord"][i], cell["magmom"][i]) for i in 1:cell["size"]]
        else
            sites = [(cell["atom"][i], cell["frac_coord"][i], cell["magmom"][i]) for i in 1:cell["size"]]
        end
    else
        if cartesian
            sites = [(cell["atom"][i], cell["cartesian_coord"][i]) for i in 1:cell["size"]]
        else
            sites = [(cell["atom"][i], cell["frac_coord"][i]) for i in 1:cell["size"]]
        end
    end
    return sites
end

"""
get_full_magnetic_sites(cell; cartesian::Bool=false, mag::Bool=cell["mag"])
"""
function get_full_magnetic_sites(cell; cartesian::Bool=false, mag::Bool=cell["mag"], tor=1e-5)
    if mag
        if cartesian
            sites = [(cell["atom"][i], cell["cartesian_coord"][i], cell["magmom"][i]) for i in 1:cell["size"] if !all(abs.(round.(cell["magmom"][i], digits=5)) .< tor)]
        else
            sites = [(cell["atom"][i], cell["frac_coord"][i], cell["magmom"][i]) for i in 1:cell["size"] if !all(abs.(round.(cell["magmom"][i], digits=5)) .< tor)]
        end
    else
        error("can not get full magnetic sites for non-magnetic cell")
    end
    return sites
end

"""
get_site_on_nth_cell(cell, i::Int64, n::Vector{Int64}; cartesian::Bool=false, mag::Bool=cell["mag"])
"""
function get_site_on_nth_cell(cell, i::Int64, n::Vector{Int64}; cartesian::Bool=false, mag::Bool=cell["mag"])
    lattice = cell["lattice"]
    site0 = get_a_site(cell, i;cartesian=cartesian, mag=mag)
    if mag
        if cartesian
            siten = (site0[1], site0[2]+lattice*n, site0[3])
        else
            siten = (site0[1], site0[2]+n, site0[3])
        end
    else
        if cartesian
            siten = (site0[1], site0[2]+lattice*n)
        else
            siten = (site0[1], site0[2]+n)
        end
    end
end

"""
get_index_and_cell_of_a_site(cell, site; cartesian::Bool=false) -> Tuple{Int64, Vector{Int64}}
"""
function get_index_and_cell_of_a_site(cell, site; cartesian::Bool=false)::Tuple{Int64, Vector{Int64}}
    for i in 1:cell["size"]
        sitei = get_a_site(cell, i, cartesian=cartesian)
        if is_site_same(site, sitei)
            if cartesian
                error()
            else
                return i, Int64.(round.(site[2]-sitei[2], digits=3))
            end
        end
    end
    error()
end

"""
get_a_site(cell, i::Int64; cartesian::Bool=false, mag::Bool=cell["mag"])
"""
function get_a_site(cell, i::Int64; cartesian::Bool=false, mag::Bool=cell["mag"])
    if mag
        if cartesian
            site = (cell["atom"][i], cell["cartesian_coord"][i], cell["magmom"][i])
        else
            site = (cell["atom"][i], cell["frac_coord"][i], cell["magmom"][i])
        end
    else
        if cartesian
            site = (cell["atom"][i], cell["cartesian_coord"][i])
        else
            site = (cell["atom"][i], cell["frac_coord"][i])
        end
    end
    return site
end

"""
is_site_same(i, j, cell; mag::Bool) -> Bool
is_site_same(site1, site2) -> Bool
"""
function is_site_same(cell, i::Int64, j::Int64; mag::Bool=cell["mag"])::Bool
    sitei = get_a_site(cell, i, mag=mag)
    sitej = get_a_site(cell, j, mag=mag)
    return is_site_same(sitei, sitej)
end
function is_site_same(site1, site2)::Bool
    if length(site1) == length(site2) == 3
        (atom1, coord1, magmom1) = site1
        (atom2, coord2, magmom2) = site2
        if is_atom_same(atom1, atom2) && is_coord_same(coord1, coord2) && is_magmom_same(magmom1, magmom2)
            return true
        else
            return false
        end
    elseif length(site1) == length(site2) == 2
        (atom1, coord1) = site1
        (atom2, coord2) = site2
        if is_atom_same(atom1, atom2) && is_coord_same(coord1, coord2)
            return true
        else
            return false
        end
    end
    error()
end

"""
is_atom_same(atom1, atom2) -> Bool
"""
function is_atom_same(atom1, atom2)::Bool
    if atom1 == atom2
        return true
    else
        return false
    end
end

"""
is_coord_same(coord1, coord2; tor=1e-3) -> Bool
"""
function is_coord_same(coord1, coord2; tor=1e-3)::Bool
    relative_coord = mod.(round.(coord2-coord1, digits=2), 1)
    if norm(relative_coord) < tor
        #println(coord2)
        #println(coord1)
        #println()
        return true
    else
        return false
    end
end

"""
is_magmom_same(magmom1, magmom2; tor=1e-3) -> Bool
"""
function is_magmom_same(magmom1, magmom2; tor=1e-3)::Bool
    relative_magmom = round.(magmom1-magmom2, digits=5)
    if norm(relative_magmom) < tor
        return true
    else
        return false
    end
end

end#(module SystemSite)

###################################################SystemSymmetry#####################################################
######################################################################################################################
module SystemSymmetry
using LinearAlgebra
using ..SystemStructure
using ..SystemSite
using ...RealSpace
using ...Symmetry
#using ...SpinPointGroup:get_class_of_spin_point_group, spg_symtype, spg_data
#using ...PointGroup:get_isometry_type

export get_symmetries
export get_full_spin_space_group_symmetries
export get_full_spin_space_group_symmetries!
export get_full_spin_point_group_symmetries
export get_full_spin_point_group_symmetries!
export get_magnetic_space_group_symmetries
export get_magnetic_space_group_symmetries!
export get_magnetic_point_group_symmetries
export get_magnetic_point_group_symmetries!
export get_spin_rotation
export get_point_group_symmetries
export get_point_group_symmetries!
export get_space_group_symmetries
export get_space_group_symmetries!
export get_translation_group_symmetries
export get_translation
export get_lattice_symmetries
export get_lattice_symmetries!
export get_metric_of_lattice
export isa_symmetry
export perform_operation_on_site

"""
get_symmetries(cell)
"""
function get_symmetries(cell)
    if cell["mag"]
        if cell["soc"]
            symmetry = get_magnetic_space_group_symmetries(cell)
        else
            symmetry = get_full_spin_space_group_symmetries(cell)
        end
    else
        symmetry = get_space_group_symmetries(cell)
    end
    return symmetry
end

"""
get_full_spin_space_group_symmetries(cell) -> Dict{String, Any}
get_full_spin_space_group_symmetries(cell)! -> Dict{String, Any}
"""
function get_full_spin_space_group_symmetries(cell)::Dict{String, Any}
    if cell["magmom_arrangement"] in [3, 4]
        println("Attention: may not find all symmetries for coplanar or noncoplanar magmom arrangement due to the infinity of magnetic translation groups")
    end
    mats = Vector{spinspasym}()
    for (r, t) in get_space_group_symmetries(cell)["operations"]
        for s in get_spin_rotation(cell, (r=r, t=t), soc=false)
            push!(mats, (r=r, t=t, s=s))
        end
    end
    s0 = get_one_pure_spin_rotation(cell)
    if s0 != nothing
        mat0 = (r=Matrix{Int64}(I,3,3), t=[0.0, 0.0, 0.0], s=s0)
        mats = append!(mats, [composite_symmetries(mat0, mat, modtrans=false) for mat in mats])
    end
    if !isa_group(mats)
        for mat in mats println(mat) end
        error("Error occurs at get_full_spin_space_group_symmetries: symmetries don't form a group")
    end
    spinspace_syms = Dict("size"=>length(mats), "operations"=>mats)
    return spinspace_syms
end
function get_full_spin_space_group_symmetries!(cell)::Dict{String, Any}
    spinspace_syms = get_full_spin_space_group_symmetries(cell)
    cell["ssg"] = spinspace_syms
    return spinspace_syms
end

#=
"""
get_full_spin_point_group_symmetries(cell) -> Dict{String, Any}
get_full_spin_point_group_symmetries(cell)! -> Dict{String, Any}
"""
function get_full_spin_point_group_symmetries(cell)::Dict{String, Any}
    if cell["magmom_arrangement"] in [3, 4]
        println("Attention: may not find all symmetries for coplanar or noncoplanar magmom arrangement due to the choice of spin coordination")
    end
    mats = Vector{spinspasym}()
    for (r, t) in get_point_group_symmetries(cell)["operations"]
        for s in get_spin_rotation(cell, (r, t), soc=false)
            push!(mats, (r=r, t=t, s=s))
        end
    end
    s0 = get_one_pure_spin_rotation(cell)
    if s0 != nothing
        mat0 = (r=Matrix{Int64}(I,3,3), t=[0.0, 0.0, 0.0], s=s0)
        mats = append!(mats, [composite_symmetries(mat0, mat, modtrans=false) for mat in mats])
    end
    if !isa_group(mats) error("Error occurs at get_full_spin_point_group_symmetries: symmetries don't form a group") end
    #spinpoint_syms = Dict("size"=>length(mats), "operations"=>mats)
    nsyms = Vector{spg_symtype}()
    for mat in mats
       push!(nsyms, (r=mat.r, s=Int64.(mat.s)))
    end
    #fspgidx = get_class_of_spin_point_group(syms, spintype = cell["magmom_arrangement"])
    fspgidx = get_class_of_spin_point_group(nsyms, SO2=3)
    nspgidx = spg_data[fspgidx].nspg
    spinpoint_syms = Dict("size"=>length(mats), "operations"=>mats, "fspgindex"=>fspgidx, "nspgindex"=>nspgidx)
    return spinpoint_syms
end
function get_full_spin_point_group_symmetries!(cell)::Dict{String, Any}
    spinpoint_syms = get_full_spin_point_group_symmetries(cell)
    cell["fspg"] = spinpoint_syms["fspgindex"]
    cell["nspg"] = spinpoint_syms["nspgindex"]
    return spinpoint_syms
end
=#

"""
get_magnetic_space_group_symmetries(cell) -> Dict{String, Any}
get_magnetic_space_group_symmetries!(cell) -> Dict{String, Any}
"""
function get_magnetic_space_group_symmetries(cell)::Dict{String, Any}
    syms = Vector{spinspasym}()
    for (r, t) in get_space_group_symmetries(cell)["operations"]
        for s in get_spin_rotation(cell, (r=r, t=t), soc=true)
            #T = Int64(round(det(s), digits=5))
            push!(syms, (r=r, t=t, s=s))
        end
    end
    if !isa_group(syms) error("Error occurs at get_magnetic_space_group_symmetries: symmetries don't form a group") end
    magspace_syms = Dict("size"=>length(syms), "operations"=>syms)
    return magspace_syms
end
function get_magnetic_space_group_symmetries!(cell)::Dict{String, Any}
    magspace_syms = get_magnetic_space_group_symmetries(cell)
    cell["msg"] = magspace_syms
    return magspace_syms
end

"""
get_magnetic_point_group_symmetries(cell) -> Dict{String, Any}
get_magnetic_point_group_symmetries(cell)! -> Dict{String, Any}
"""
function get_magnetic_point_group_symmetries(cell)::Dict{String, Any}
    syms = Vector{spinspasym}()
    for (r, t) in get_point_group_symmetries(cell)["operations"]
        for s in get_spin_rotation(cell, (r, t), soc=true)
            push!(syms, (r=r, t=t, s=s))
        end
    end
    if !isa_group(syms) error("Error occurs at get_magnetic_point_group_symmetries: symmetries don't form a group") end
    magpoint_syms = Dict("size"=>length(syms), "operations"=>syms)
    return magpoint_syms
end
function get_magnetic_point_group_symmetries!(cell)::Dict{String, Any}
    magpoint_syms = get_magnetic_point_group_symmetries(cell)
    cell["mpg"] = magpoint_syms
    return magpoint_syms
end

"""
get_spin_rotation(cell, (r, t)::spasym; soc::Bool) -> Vector{Matrix{Float64}}
get_spin_rotation(cell, r::ptsym; soc::Bool) -> Vector{Matrix{Float64}}
"""
function get_spin_rotation(cell, (r, t); soc::Bool)::Vector{Matrix{Float64}}
    #=
    Warning: only get spin rotation belonging to 32 point group
    =#
    relative_axes = [
        [ 1,  0,  0], [ 0,  1,  0], [ 0,  0,  1], [-1,  0,  0], [ 0, -1,  0],
        [ 0,  0, -1], [ 0,  1,  1], [ 1,  0,  1], [ 1,  1,  0], [ 0, -1, -1],
        [-1,  0, -1], [-1, -1,  0], [ 0,  1, -1], [-1,  0,  1], [ 1, -1,  0],
        [ 0, -1,  1], [ 1,  0, -1], [-1,  1,  0], [ 1,  1,  1], [-1, -1, -1],
        [-1,  1,  1], [ 1, -1,  1], [ 1,  1, -1], [ 1, -1, -1], [-1,  1, -1],
        [-1, -1,  1]
    ]
    magtype = cell["magmom_arrangement"]
    spin_rotations = Vector{Matrix{Float64}}()
    if soc
        for mat in [r, -r]
            mat = Float64.(mat)
            sym = (r=r, t=t, s=mat)
            if isa_symmetry(sym, cell, mag=true) push!(spin_rotations, mat) end
        end
    else
        if magtype == 0
            mat = Matrix{Float64}(I,3,3)
            push!(spin_rotations, mat)
        elseif magtype == 1
            for mat in [Matrix{Float64}(I,3,3)]
                sym = (r=r, t=t, s=mat)
                if isa_symmetry(sym, cell, mag=true) push!(spin_rotations, mat) end
            end
        elseif magtype == 2
            for mat in [Matrix{Float64}(I,3,3), -Matrix{Float64}(I,3,3)]
                sym = (r=r, t=t, s=mat)
                if isa_symmetry(sym, cell, mag=true) push!(spin_rotations, mat) end
            end
        elseif magtype == 3
            for i in 1:26, j in 1:26, k in 1:26
                mat = [relative_axes[i];;relative_axes[j];;relative_axes[k]]
                #if !(det(mat) in [1, -1]) || get_notation_of_symmetry(mat) == nothing continue end
                if det(mat) != 1 || get_notation_of_symmetry(mat) == nothing continue end
                order = get_symmetry_order_of_point_symmetry(mat)
                if composite_symmetries(mat, order=order) != Matrix{Int64}(I,3,3) ||
                    !(transpose(mat)*mat in [
                    [1  0 0;  0 1 0; 0 0 1],
                    [1 -1 0; -1 2 0; 0 0 1],
                    [2 -1 0; -1 1 0; 0 0 1]
                ])
                    continue
                end
                mat = Float64.(mat)
                sym = (r=r, t=t, s=mat)
                if isa_symmetry(sym, cell, mag=true) push!(spin_rotations, mat) end
            end
        elseif magtype == 4
            #=warning: still under modified=#
            for i in 1:26, j in 1:26, k in 1:26
                mat = [relative_axes[i];;relative_axes[j];;relative_axes[k]]
                if !(det(mat) in [1, -1]) continue end
                mat = Float64.(mat)
                sym = (r=r, t=t, s=mat)
                if isa_symmetry(sym, cell, mag=true) push!(spin_rotations, mat) end
            end
            if length(spin_rotations) >= 2 error("Error occurs at get_spin_rotation: noncoplanar can not have pure spin rotation") end
        end
    end
    return spin_rotations
end

"""
get_one_pure_spin_rotation(cell)::Union{Nothing, Matrix{Float64}}
"""
function get_one_pure_spin_rotation(cell)#::Union{Nothing, Matrix{Float64}}
    magtype = cell["magmom_arrangement"]
    if !(magtype in [1, 2, 3]) return nothing end
    r = Matrix{Int64}(I,3,3)
    t = [0.0, 0.0, 0.0]
    for mat in [diagm([-1.0, 1.0, 1.0]), diagm([1.0, -1.0, 1.0]), diagm([1.0, 1.0, -1.0])]
        sym = (r=r, t=t, s=mat)
        if isa_symmetry(sym, cell, mag=true) return mat end
    end
    error("can not get one pure spin rotations")
end

"""
get_point_group_symmetries(cell) -> Dict{String, Any}
get_point_group_symmetries!(cell) -> Dict{String, Any}
"""
function get_point_group_symmetries(cell; tor=1e-3)#::Dict{String, Any}
    symmorphic_space_syms = Dict(
        0 => Vector{Tuple{spasym, rtype}}(),
        1 => Vector{Tuple{spasym, rtype}}(),
        2 => Vector{Tuple{spasym, rtype}}(),
        3 => Vector{Tuple{spasym, rtype}}()
    )

    for sym in get_space_group_symmetries(cell)["operations"]
        if isa_pointsymmetry(sym)
            fixedpoints = get_fixedpoints_of_symmetry_operation(sym)
            geotype = get_r_geometric_type(fixedpoints[1][2])
            for (nsym, fixedpoint) in fixedpoints
                push!(symmorphic_space_syms[geotype], (nsym, fixedpoint))
            end
        end
    end

    point_syms_list = Vector{Dict{String, Any}}()
    for (sym, fixedpoint) in symmorphic_space_syms[0]
        for i in eachindex(point_syms_list)
            ptsyms = point_syms_list[i]
            if is_r_same(fixedpoint, ptsyms["fixedpoint"], modtrans=false)
                push!(point_syms_list[i]["operations"], sym)
                @goto next_sym0
            end
        end
        push!(point_syms_list, Dict("geotype"=>0, "fixedpoint"=>fixedpoint, "operations"=>[sym]))
        @label next_sym0
    end

    for (sym, fixedpoint) in symmorphic_space_syms[1]
        stat = false
        for i in eachindex(point_syms_list)
            ptsyms = point_syms_list[i]
            if ptsyms["geotype"] == 0
                if is_r1_contained_in_r2(r1=ptsyms["fixedpoint"], r2=fixedpoint)
                    stat = true
                    push!(point_syms_list[i]["operations"], sym)
                end
            elseif ptsyms["geotype"] == 1
                intersec = get_intersection_of_r(ptsyms["fixedpoint"], fixedpoint)
                if intersec != nothing
                    geotype = get_r_geometric_type(intersec)
                    point_syms_list[i]["geotype"] = geotype
                    point_syms_list[i]["fixedpoint"] = intersec
                    stat = true
                    push!(point_syms_list[i]["operations"], sym)
                end
            end
        end
        if !stat push!(point_syms_list, Dict("geotype"=>1, "fixedpoint"=>fixedpoint, "operations"=>[sym])) end
    end

    for (sym, fixedpoint) in symmorphic_space_syms[2]
        stat = false
        for i in eachindex(point_syms_list)
            ptsyms = point_syms_list[i]
            if ptsyms["geotype"] in [0, 1]
                if is_r1_contained_in_r2(r1=ptsyms["fixedpoint"], r2=fixedpoint)
                    stat = true
                    push!(point_syms_list[i]["operations"], sym)
                end
            elseif ptsyms["geotype"] == 2
                intersec = get_intersection_of_r(ptsyms["fixedpoint"], fixedpoint)
                if intersec != nothing
                    geotype = get_r_geometric_type(intersec)
                    point_syms_list[i]["geotype"] = geotype
                    point_syms_list[i]["fixedpoint"] = intersec
                    stat = true
                    push!(point_syms_list[i]["operations"], sym)
                end
            end
        end
        if !stat push!(point_syms_list, Dict("geotype"=>2, "fixedpoint"=>fixedpoint, "operations"=>[sym])) end
    end

    for (sym, r) in symmorphic_space_syms[3]
        for i in eachindex(point_syms_list)
            push!(point_syms_list[i]["operations"], sym)
        end
    end

    point_syms = Dict("size"=>0, "fixedpoint"=>[], "operations"=>[])
    for i in eachindex(point_syms_list)
        syms = point_syms_list[i]["operations"]
        size = length(point_syms_list[i]["operations"])
        if point_syms["size"] < size
            point_syms["size"] = size
            point_syms["fixedpoint"] = point_syms_list[i]["fixedpoint"]
            point_syms["operations"] = point_syms_list[i]["operations"]
        elseif point_syms["size"] == size
            fixed1 = point_syms["fixedpoint"]
            fixed2 = point_syms_list[i]["fixedpoint"]
            for j in 1:3
                if abs(round(fixed2[j, 4]-fixed1[j, 4], digits=3)) < tor
                    continue
                else
                    if abs.(round.(abs(fixed2[j, 4])-abs(fixed1[j, 4]), digits=3)) < tor
                        if fixed2[j, 4] > fixed1[j, 4]
                            point_syms["size"] = size
                            point_syms["fixedpoint"] = point_syms_list[i]["fixedpoint"]
                            point_syms["operations"] = point_syms_list[i]["operations"]
                        end
                    elseif round.(abs(fixed2[j, 4])-abs(fixed1[j, 4]), digits=3) < tor
                        point_syms["size"] = size
                        point_syms["fixedpoint"] = point_syms_list[i]["fixedpoint"]
                        point_syms["operations"] = point_syms_list[i]["operations"]
                    end
                    break
                end
            end
        end
    end
    syms = point_syms["operations"]
    if !isa_group(syms) error("Error occurs at get_point_group_symmetries: symmetries don't form a group") end

    return point_syms
end
function get_point_group_symmetries!(cell)::Dict{String, Any}
    point_syms = get_point_group_symmetries(cell)
    cell["pg"] = point_syms
    return point_syms
end

"""
get_space_group_symmetries(cell) -> Dict{String, Any}
get_space_group_symmetries!(cell) -> Dict{String, Any}
#not considering magnetic moments of cell, get full space group symmetries
"""
function get_space_group_symmetries(cell)::Dict{String, Any}
    syms = Vector{spasym}()
    if haskey(cell, "latticesymmetry")
        lat_syms = cell["latticesymmetry"]["operations"]
    else
        lat_syms = get_lattice_symmetries!(cell)["operations"]
    end
    for r in lat_syms
        for t in get_translation(cell, r)
            push!(syms, (r=r, t=round.(t, digits=5)))
        end
    end
    if !isa_group(syms)
        for sym in syms
            println(sym)
        end
        error("Error occurs at get_space_group_symmetries: symmetries don't form a group")
    end
    space_syms = Dict("size"=>length(syms), "operations"=>syms)
    return space_syms
end
function get_space_group_symmetries!(cell)::Dict{String, Any}
    space_syms = get_space_group_symmetries(cell)
    cell["sg"] = space_syms
    return space_syms
end


"""
get_translation_group_symmetries(cell) -> Dict{String, Any}
"""
function get_translation_group_symmetries(cell)::Dict{String, Any}
    syms = Vector{spasym}()
    Isym  = [1 0 0;0 1 0;0 0 1]
    for t in get_translation(cell, Isym)
        push!(syms, (r=Isym, t=t))
    end
    if !isa_group(syms) error("Error occurs at get_translation_group_symmetries: symmetries don't form a group") end
    trans_syms = Dict("size"=>length(syms), "operations"=>syms)
    return trans_syms
end

"""
get_translation(cell, r::ptsym) -> Vector{Vector{Float64}}
#not considering magnetic moments of cell, get full translations corresonding to a point symmetry
"""
function get_translation(cell, r::ptsym; tor=1e-3)::Vector{Vector{Float64}}
    translations = Vector{Vector{Float64}}()
    coords = cell["frac_coord"]
    coord0 = r*coords[1]
    for coord in coords
        t = mod.(round.(coord-coord0, digits=4), 1)
        sym = (r=r, t=t)
        if !isa_symmetry(sym, cell, mag=false) continue end
        push!(translations, t)
        if isa_pointsymmetry(sym) continue end
        for x in [-2.0, -1.0, 0.0, 1.0, 2.0], y in [-2.0, -1.0, 0.0, 1.0, 2.0], z in [-2.0, -1.0, 0.0, 1.0, 2.0]
            sym = (r=r, t=t+[x, y, z])
            if isa_pointsymmetry(sym)
                translations[end] = t+[x, y, z]
                break
            end
        end
    end
    return translations
end

"""
get_lattice_symmetries(cell) -> Dict{String, Any}
get_lattice_symmetries!(cell) -> Dict{String, Any}
"""
function get_lattice_symmetries(cell; tor=1e-3)::Dict{String, Any}
    relative_axes = [
        [ 1,  0,  0], [ 0,  1,  0], [ 0,  0,  1], [-1,  0,  0], [ 0, -1,  0],
        [ 0,  0, -1], [ 0,  1,  1], [ 1,  0,  1], [ 1,  1,  0], [ 0, -1, -1],
        [-1,  0, -1], [-1, -1,  0], [ 0,  1, -1], [-1,  0,  1], [ 1, -1,  0],
        [ 0, -1,  1], [ 1,  0, -1], [-1,  1,  0], [ 1,  1,  1], [-1, -1, -1],
        [-1,  1,  1], [ 1, -1,  1], [ 1,  1, -1], [ 1, -1, -1], [-1,  1, -1],
        [-1, -1,  1]
    ]
    lat_ori = cell["lattice"]
    syms = Vector{ptsym}()
    metric_ori = get_metric_of_lattice(lat_ori)
    for i in 1:26, j in 1:26, k in 1:26
        mat = [relative_axes[i];;relative_axes[j];;relative_axes[k]]
        if !(det(mat) in [1, -1]) continue end
        lat = lat_ori*mat
        metric = get_metric_of_lattice(lat)
        if all(abs.(round.(metric_ori-metric, digits=3)) .< tor)
            push!(syms,mat)
        end
    end
    if !isa_group(syms) error("Error occurs at get_lattice_symmetries: symmetries don't form a group") end
    lat_syms = Dict("size"=>length(syms), "operations"=>syms)
    return lat_syms
end
function get_lattice_symmetries!(cell)::Dict{String, Any}
    lat_syms = get_lattice_symmetries(cell)
    cell["latticesymmetry"] = lat_syms
    return lat_syms
end

"""
get_metric_of_lattice(lat) -> Matrix{Float64}
"""
function get_metric_of_lattice(lat)::Matrix{Float64}
    metric = transpose(lat)*lat
    return metric
end

"""
isa_symmetry(sym::T, cell; mag::Bool) -> Bool where {T<:Union{ptsym, spasym, spinspasym}}
"""
function isa_symmetry(sym::T, cell; mag::Bool)::Bool where {T<:Union{ptsym, spasym, spinspasym}}
    sites = get_full_sites(cell, mag=mag)
    for site_ori in sites
        site_rot = perform_operation_on_site(sym, site_ori)
        for site in sites
            if is_site_same(site_rot, site) @goto next end
        end
        return false
        @label next
    end
    return true
end

end#(module SystemSymmetry)

###################################################SystemNeighbor#####################################################
######################################################################################################################
module SystemNeighbor
using LinearAlgebra
using ...MathTool
using ..SystemSite
using ..SystemSymmetry

export get_neighbors
export get_neighbors_between_atoms
export get_neighbors_of_certain_distance
export get_neighbors_of_certain_distance_between_atoms
export get_neighbor_distances
export get_neighbor_distances_between_atoms
export sort_neighbors_into_equivalent_ones
export get_equivalent_neighbor
export is_neighbor_contained

#=
get_neighbors(cell, n::int64) -> Dict{Int64, Any}
get_neighbors(cell, range::UnitRange{Int64}) -> Dict{Int64, Any}
get_neighbors_of_certain_distance(cell, dist0::Float64; tor=1e-3)
get_neighbor_distances(cell, n::Int64; tor=1e-3) -> Vector{Float64}
get_neighbor_distances(cell, range::UnitRange{Int64}; tor=1e-3) -> Vector{Float64}
get_neighbors_between_atoms(cell, n::Int64, atom1, atom2) -> Dict{Int64, Any}
get_neighbors_between_atoms(cell, range::UnitRange{Int64}, atom1, atom2) -> Dict{Int64, Any}
get_neighbors_of_certain_distance_between_atoms(cell, dist0::Float64, atom1, atom2; tor=1e-3)
get_neighbor_distances_between_atoms(cell, n::Int64, atom1, atom2; tor=1e-3) -> Vector{Float64}
get_neighbor_distances_between_atoms(cell, range::UnitRange{Int64}, atom1, atom2; tor=1e-3) -> Vector{Float64}
sort_neighbors_into_equivalent_ones(cell, neighbors::Vector{Any}, mats) -> Dict{Int64, Any}
get_equivalent_neighbor(cell, neighbor, mats) -> Vector{Any}
is_neighbor_contained(neighbor0, neighbors) -> Bool
=#


"""
get_neighbors(cell, n::int64) -> Dict{Int64, Any}
get_neighbors(cell, range::UnitRange{Int64}) -> Dict{Int64, Any}
"""
function get_neighbors(cell, n::Int64; sym::Bool=false)::Dict{Int64, Dict{Any, Any}}
    if n < 1 error("Error occurs at get_neighbors:please enter positive number n") end
    return get_neighbors(cell, 1:n, sym=sym)
end
function get_neighbors(cell, range::UnitRange{Int64}; sym::Bool=false)::Dict{Int64, Dict{Any, Any}}
    dists = get_neighbor_distances(cell, range)
    #neighbors = Dict{Int64, Any}()
    neighbors = Dict{Int64, Dict{Any, Any}}()
    for i in eachindex(dists)
        #neighbors[i] = get_neighbors_of_certain_distance(cell, dists[i])
        neighbors[i] = Dict("dist" => dists[i], 1 => get_neighbors_of_certain_distance_between_atoms(cell, dists[i], atom1, atom2))
    end
    if sym
        mats = get_symmetries(cell)["operations"]
        sort_neighbors_into_equivalent_ones!(cell, neighbors, mats)
    end
    return neighbors
end

"""
get_neighbors_of_certain_distance(cell, dist0::Float64; tor=1e-3)
"""
function get_neighbors_of_certain_distance(cell, dist0::Float64; tor=1e-3)
    neighbors = []
    size = cell["size"]
    bonds = [(i, j) for i=1:size for j=1:size]
    relative_cell = [[i, j, k] for i=-2:2 for j=-2:2 for k=-2:2]
    for (i, j) in bonds
        for n in relative_cell
            sitei = get_site_on_nth_cell(cell, i, [0, 0, 0], cartesian=true)
            sitej = get_site_on_nth_cell(cell, j, n, cartesian=true)
            dist = get_distance(sitei[2], sitej[2])
            if abs(round(dist-dist0, digits=4)) > tor continue end
            push!(neighbors, [i, j, n])
        end
    end
    return neighbors
end

"""
get_neighbor_distances(cell, n::Int64; tor=1e-3) -> Vector{Float64}
get_neighbor_distances(cell, range::UnitRange{Int64}; tor=1e-3) -> Vector{Float64}
"""
function get_neighbor_distances(cell, n::Int64; tor=1e-3)::Vector{Float64}
    if n < 1 error("Error occurs at get_neighbor_distances:please enter positive number n") end
    return get_neighbor_distances(cell, 1:n)
end
function get_neighbor_distances(cell, range::UnitRange{Int64}; tor=1e-3)::Vector{Float64}
    dists = Float64.([100+i for i=0:range[end]])
    lattice = cell["lattice"]
    size = cell["size"]
    bonds = [(i, j) for i=1:size for j=i:size]
    relative_cell = [[i, j, k] for i=-2:2 for j=-2:2 for k=-2:2]
    for (i, j) in bonds
        for n in relative_cell
            sitei0 = get_site_on_nth_cell(cell, i, [0, 0, 0], cartesian=true)
            sitejn = get_site_on_nth_cell(cell, j, n, cartesian=true)
            dist = get_distance(sitei0[2], sitejn[2])
            for x in dists
                if abs(round(dist-x, digits=4)) < tor @goto next end
            end
            push!(dists, dist)
            pop!(sort!(dists))
            @label next
        end
    end
    popat!(dists, 1)
    return dists[range]
end

"""
get_neighbors_between_atoms(cell, n::Int64, atom1, atom2) -> Dict{Int64, Dict{Any, Any}}
get_neighbors_between_atoms(cell, range::UnitRange{Int64}, atom1, atom2) -> Dict{Int64, Dict{Any, Any}}
"""
function get_neighbors_between_atoms(cell, n::Int64, atom1, atom2; sym::Bool=false)::Dict{Int64, Dict{Any, Any}}
    if n < 1 error("Error occurs at get_neighbors_between_atoms:please enter positive number n") end
    return get_neighbors_between_atoms(cell, 1:n, atom1, atom2, sym=sym)
end
function get_neighbors_between_atoms(cell, range::UnitRange{Int64}, atom1, atom2; sym::Bool=false)::Union{Dict{Int64, Any}, Dict{Int64, Dict{Any, Any}}}
    dists = get_neighbor_distances_between_atoms(cell, range, atom1, atom2)
    #neighbors = Dict{Int64, Any}()
    neighbors = Dict{Int64, Dict{Any, Any}}()
    for i in eachindex(dists)
        neighbors[i] = Dict("dist" => dists[i], 1 => get_neighbors_of_certain_distance_between_atoms(cell, dists[i], atom1, atom2))
        #neighbors[i] = get_neighbors_of_certain_distance_between_atoms(cell, dists[i], atom1, atom2)
    end
    if sym
        mats = get_symmetries(cell)["operations"]
        sort_neighbors_into_equivalent_ones!(cell, neighbors, mats)
    end
    return neighbors
end

"""
get_neighbors_of_certain_distance_between_atoms(cell, dist0::Float64, atom1, atom2; tor=1e-3)
Return [[i, j, n], ...]
i start atom
j end atom
n vector from cell of atomi to cell of atomj
"""
function get_neighbors_of_certain_distance_between_atoms(cell, dist0::Float64, atom1, atom2; tor=1e-3)
    neighbors = []
    size = cell["size"]
    bonds = [(i, j) for i=1:size for j=1:size]
    relative_cell = [[i, j, k] for i=-2:2 for j=-2:2 for k=-2:2]
    for (i, j) in bonds
        for n in relative_cell
            sitei = get_site_on_nth_cell(cell, i, [0, 0, 0], cartesian=true)
            sitej = get_site_on_nth_cell(cell, j, n, cartesian=true)
            dist = get_distance(sitei[2], sitej[2])
            if abs(round(dist-dist0, digits=4)) > tor continue end
            if !(sitei[1] == atom1 && sitej[1] == atom2) && !(sitei[1] == atom2 && sitej[1] == atom1) continue end
            push!(neighbors, [i, j, n])
        end
    end
    return neighbors
end

"""
get_neighbor_distances_between_atoms(cell, n::Int64, atom1, atom2; tor=1e-3) -> Vector{Float64}
get_neighbor_distances_between_atoms(cell, range::UnitRange{Int64}, atom1, atom2; tor=1e-3) -> Vector{Float64}
"""
function get_neighbor_distances_between_atoms(cell, n::Int64, atom1, atom2; tor=1e-3)::Vector{Float64}
    if n < 1 error("Error occurs at get_neighbor_distances_between_atoms:please enter positive number n") end
    return get_neighbor_distances_between_atoms(cell, 1:n, atom1, atom2)
end
function get_neighbor_distances_between_atoms(cell, range::UnitRange{Int64}, atom1, atom2; tor=1e-3)::Vector{Float64}
    dists = Float64.([100+i for i in 0:range[end]])
    lattice = cell["lattice"]
    size = cell["size"]
    bonds = [(i, j) for i=1:size for j=i:size if cell["atom"][i] == atom1 && cell["atom"][j] == atom2]
    relative_cell = [[i, j, k] for i=-2:2 for j=-2:2 for k=-2:2]
    for (i, j) in bonds
        for n in relative_cell
            sitei0 = get_site_on_nth_cell(cell, i, [0, 0, 0], cartesian=true)
            sitejn = get_site_on_nth_cell(cell, j, n, cartesian=true)
            dist = get_distance(sitei0[2], sitejn[2])
            for x in dists
                if abs(round(dist-x, digits=4)) < tor @goto next end
            end
            push!(dists, dist)
            pop!(sort!(dists))
            @label next
        end
    end
    popat!(dists, 1)
    return dists[range]
end

"""
sort_neighbors_into_equivalent_ones(cell, neighbors::Vector{Any}, mats) -> Dict{Any, Any}
"""
function sort_neighbors_into_equivalent_ones!(cell, neighbors_list::Dict{Int64, Dict{Any, Any}}, mats)::Dict{Int64, Dict{Any, Any}}
    neighbors_list_ori = deepcopy(neighbors_list)
    for (i, neighbors) in neighbors_list_ori
        neighbors_list[i] = sort_neighbors_into_equivalent_ones(cell, neighbors[1], mats)
        neighbors_list[i]["dist"] = neighbors["dist"]
        #neighbors_list[i] = sort_neighbors_into_equivalent_ones(cell, neighbors, mats)
    end
    return neighbors_list
end
function sort_neighbors_into_equivalent_ones(cell, neighbors::Vector{Any}, mats)::Dict{Any, Any}
    neighbors_ori = deepcopy(neighbors)
    neighbors = Dict{Int64, Any}()
    n = 1
    for neighbor in neighbors_ori
        neighbor_inv = [neighbor[2], neighbor[1], flipsign.(neighbor[3], -1)]
        for eq_neighbors in values(neighbors)
            if is_neighbor_contained(neighbor, eq_neighbors) || is_neighbor_contained(neighbor_inv, eq_neighbors) @goto next_neighbor end
        end
        eq_neighbors = get_equivalent_neighbor(cell, neighbor, mats)
        neighbors[n] = eq_neighbors
        n += 1
        @label next_neighbor
    end
    return neighbors
end

"""
get_equivalent_neighbor(cell, neighbor, mats) -> Vector{Any}
"""
function get_equivalent_neighbor(cell, neighbor, mats)::Vector{Any}
    eq_neighbors = [neighbor]
    i = neighbor[1]
    j = neighbor[2]
    cellj = neighbor[3]
    for mat in mats
        k, cellk = perform_operation_on_site(cell, mat, (i, [0, 0, 0]))
        l, celll = perform_operation_on_site(cell, mat, (j, cellj))
        neighbor_rot = [k, l, celll-cellk]
        neighbor_rotinv = [l, k, cellk-celll]
        #println(neighbor_rot, is_neighbor_contained(neighbor_rot, eq_neighbors))
        if !is_neighbor_contained(neighbor_rot, eq_neighbors) push!(eq_neighbors, neighbor_rot) end
        if !is_neighbor_contained(neighbor_rotinv, eq_neighbors) push!(eq_neighbors, neighbor_rotinv) end
    end
    return eq_neighbors
end

"""
is_neighbor_contained(neighbor0, neighbors) -> Bool
"""
function is_neighbor_contained(neighbor0, neighbors)::Bool
    for neighbor in neighbors
        if all(neighbor0 .== neighbor) return true end
    end
    return false
end

end#(module SystemNeighbor)

###################################################SystemSymmetry#####################################################
######################################################################################################################
module kSpaceSymmetry
using LinearAlgebra
using ...kSpace
using ...Symmetry

export get_little_group_of_k

"""
get_little_group_of_k(syms; k)
"""
function get_little_group_of_k(syms; k)
    littlegroup = []
    for rsym in syms
        ksym = get_symmetry_in_kspace(rsym)
        if is_k_invariant_under_symmetry_operation(k, ksym)
            push!(littlegroup, rsym)
        end
    end
    return littlegroup
end

end#module kSpaceSymmetry

###################################################GroupRepresentation#####################################################
######################################################################################################################
module SystemGroupRepresentation
using LinearAlgebra
using ..SystemSymmetry
using ..SystemSite
using ..kSpaceSymmetry
using ...kSpace
using ...Symmetry#: spinspasym, spinptsym, spasym, ptsym, get_symmetry_in_kspace
using ...GroupIrrepresentationConstruct:get_coirreps_of_little_group_of_k
using ...GroupRepresentation:get_multiplicity_of_irrep

export get_irrep_decomposition_at_k
export get_k_space_rep_of_symmetry
export get_real_space_rep_of_symmetry


function get_irrep_decomposition_at_k(cell, mats; k)
    mats = get_little_group_of_k(mats; k=k)
    println("Constructing co-irrepresentation...")
    flush(stdout)
    Umats, irreps = get_coirreps_of_little_group_of_k(mats; k=k)
    println("Constructing representation in magnon basis...")
    flush(stdout)
    reps = get_k_space_rep_of_symmetry(cell, Umats, k=k)
    println("Getting irrep decomposition...")
    flush(stdout)
    muls = get_multiplicity_of_irrep(reps, irreps)
    println("Done.")
    for i in eachindex(muls)
        dim = Int64(size(irreps[i,1],1))
        mul = Int64(muls[i]/2)
        println(rpad("irrep$i", 8),"dim=$dim","   $mul")
    end
end

"""
get_k_space_rep_of_symmetry(cell, rsym::spinspasym; rot::Matrix{Float64}=Matrix{Float64}(I, 3, 3), k, tor=1e-5)
get_k_space_rep_of_symmetry(cell, rsym::Union{spinspasym, spinptsym}; k)
formula
    Ug
    VG*Ug*exp[-i(k+G)*t]
Return
    VG*Ug*exp[-iG*t]
"""
function get_k_space_rep_of_symmetry(cell, rsyms::Union{Vector{magspasym}, Vector{spinspasym}}; k, tor=1e-5)::Vector{Matrix{ComplexF64}}
    reps = Matrix{ComplexF64}[]
    for rsym in rsyms
        rep = get_k_space_rep_of_symmetry(cell, rsym, k=k)
        push!(reps, round.(rep, digits=3))
    end
    return reps
end
function get_k_space_rep_of_symmetry(cell, rsym::magspasym; k, tor=1e-5)::Matrix{ComplexF64}
    ksym = get_symmetry_in_kspace(rsym)
    #if !is_k_invariant_under_symmetry_operation(k, ksym) error("Error occurs at get_k_space_rep_of_symmetry: k not invariant under symmetry") end
    G = ksym * k - k
    if !all(abs.(round.(G[1:3, 1:3], digits=5)) .< tor) error("Error occurs at get_k_space_rep_of_symmetry(: can not get G=nk-k") end
    G = G[:, 4]
    Ug = get_real_space_rep_of_symmetry(cell, rsym, holespace=false)
    phase = exp(-2.0*π*im*sum(G.*rsym.t))
    sites = get_full_magnetic_sites(cell)
    VG = diagm([exp(2.0*π*im*sum(G.*sites[i][2])) for j in 1:2 for i in eachindex(sites)])
    if Int64.(rsym.s) in [[ 1 0; 0 -1], [-1 0; 0 1], [ 0 -1; 1 0], [ 0 1;-1 0]]
        Ug = [[Ug;; zeros(size(Ug))]; [zeros(size(Ug));; -conj.(Ug)]]
    else
        Ug = [[Ug;; zeros(size(Ug))]; [zeros(size(Ug));; conj.(Ug)]]
    end
    return phase * VG * Ug
end
function get_k_space_rep_of_symmetry(cell, rsym::spinspasym; k, tor=1e-5)::Matrix{ComplexF64}
    ksym = get_symmetry_in_kspace(rsym)
    #if !is_k_invariant_under_symmetry_operation(k, ksym) error("Error occurs at get_k_space_rep_of_symmetry: k not invariant under symmetry") end
    G = ksym * k - k
    if !all(abs.(round.(G[1:3, 1:3], digits=5)) .< tor)
        println(G)
        error("Error occurs at get_k_space_rep_of_symmetry: can not get G=nk-k")
    end
    G = G[:, 4]
    phase = exp(-2.0*π*im*sum(G.*rsym.t))
    sites = get_full_magnetic_sites(cell)
    VG = diagm([exp(2.0*π*im*sum(G.*sites[i][2])) for j in 1:2 for i in eachindex(sites)])
    Ug = get_real_space_rep_of_symmetry(cell, rsym)
    Ug = [[Ug;; zeros(size(Ug))]; [zeros(size(Ug));; conj.(Ug)]]
    return phase * VG * Ug
end
function get_k_space_rep_of_symmetry(cell, rsym::spinptsym; k, tor=1e-5)::Matrix{ComplexF64}
    ksym = get_symmetry_in_kspace(rsym)
    #if !is_k_invariant_under_symmetry_operation(k, ksym) error("Error occurs at get_k_space_rep_of_symmetry: k not invariant under symmetry") end
    G = ksym * k - k
    if !all(abs.(round.(G[1:3, 1:3], digits=5)) .< tor) error("Error occurs at get_k_space_rep_of_symmetry(: can not get G=nk-k") end
    G = G[1:3, 4]
    sites = get_full_magnetic_sites(cell)#get_full_sites(cell)
    VG = diagm([exp(2.0*π*im*sum(G.*sites[i][2])) for j in 1:2 for i in eachindex(sites)])
    Ug = get_real_space_rep_of_symmetry(cell, rsym)
    Ug = [[Ug;; zeros(size(Ug))]; [zeros(size(Ug));; conj.(Ug)]]
    return Vg * Ug
end

"""
get_real_space_rep_of_symmetry
rot: transformation matrix from lattice vector notation to cartesian vector notation
formula
    g a^+_{i,m} g^+ = sum_{n} a^+_{j,n} [U_g]_{nm}
    i,j : different cell
    m,n : inequivalent site in a cell
    Ug : one block of real space representation
procedure
    1.get symmetry matrix and magnetic moments in cartesian coord
    2.get rotation matrix from local_axis_to_global_axis
    3.get Ug
"""
function get_real_space_rep_of_symmetry(cell, sym::magspasym; holespace::Bool, tor=1e-3)
    sites = get_full_magnetic_sites(cell)
    cartesian_magmoms = [site[3] for site in sites]
    Rs = get_rotation_matrix_from_local_axis_to_global_axis(cartesian_magmoms)

    spinz = Rs[1]*[0,0,1]
    spinzvec = spinz/norm(spinz)
    spins = []
    for site in sites
        spin = site[3]
        spinvec = spin/norm(spin)
        cdot = Int64(round(sum(spinzvec.*spinvec), digits=5))
        push!(spins, cdot)
    end

    Ug = zeros(Int64, length(sites), length(sites))
    for i in eachindex(sites)
        sitei = sites[i]
        nsite = (sitei[1], sym.r*sitei[2]+sym.t)
        for j in eachindex(sites)
            sitej = sites[j]
            if nsite[1] == sitej[1] && all(mod.(abs.(round.(mod.(nsite[2]-sitej[2], 1), digits=4)), 1) .< tor)
                if spins[i] == 1 && spins[j] == 1
                    if !holespace
                        Ug[j, i] = sym.s[1, 1]
                    else
                        Ug[j, i] = sym.s[2, 2]
                    end
                elseif spins[i] == 1 && spins[j] == -1
                    if !holespace
                        Ug[j, i] = sym.s[1, 2]
                    else
                        Ug[j, i] = sym.s[2, 1]
                    end
                elseif spins[i] == -1 && spins[j] == 1
                    if !holespace
                        Ug[j, i] = sym.s[2, 1]
                    else
                        Ug[j, i] = sym.s[1, 2]
                    end
                elseif spins[i] == -1 && spins[j] == -1
                    if !holespace
                        Ug[j, i] = sym.s[2, 2]
                    else
                        Ug[j, i] = sym.s[1, 1]
                    end
                else
                    error("please check spin orientation")
                end
                @goto next_i
            end
        end
        error("can not find real space representation")
        @label next_i
    end
    return Ug
end
function get_real_space_rep_of_symmetry(cell, sym::Union{spinspasym, spinptsym}; tor=1e-3)
    sites = get_full_magnetic_sites(cell)
    rot = Matrix{Float64}(I, 3, 3)
    #rot = [1 -0.5 0.0;0.0 sqrt(3)/2 0.0;0.0 0.0 1.0]
    cartesian_sym = (r=sym.r, s=rot*sym.s*inv(rot))
    cartesian_magmoms = [rot*site[3] for site in sites]
    Rs = get_rotation_matrix_from_local_axis_to_global_axis(cartesian_magmoms)
    Ug = zeros(Complex{Float64}, length(sites), length(sites))
    for i in eachindex(sites)
        sitei = sites[i]
        nsite = perform_operation_on_site(sym, sitei)
        T = Int64(round(det(sym.s), digits=5))
        for j in eachindex(sites)
            sitej = sites[j]
            if nsite[1] == sitej[1] && all(mod.(abs.(round.(nsite[2]-sitej[2], digits=4)), 1) .< tor)
                factor = get_phase_factor(cartesian_sym, Rs[i], Rs[j], T)
                if factor == nothing
                    println(i," ",j," ",T)
                    println(cartesian_sym.s)
                    println(inv(Rs[i]))
                    println(inv(Rs[j]))
                    error()
                end
                Ug[j, i] = factor
                @goto next_i
            end
        end
        error("can not find real space representation")
        @label next_i
    end
    return Ug
end

"""
get_phase_factor(sym, Ri, Rj, T; tor=1e-5)
    Ri, Rj : rotation matrix from local axis to global axis of ith/jth site
formula
    operator from to local axis to global axis
    (a^x_i, a^y_i, a^z_i)^T     ->     (b^x_i, b^y_i, b^z_i)^T = Ri(a^x_i, a^y_i, a^z_i)^T
    a^+_i = a^x_i + im*a^y_i    ->     b^+_i = b^x_i + im*b^y_i

    operator under symmetry operation in global axis
    (b^x_i, b^y_i, b^z_i)^T     ->     (nb^x_i, nb^y_i, nb^z_i)^T = sym.s*(b^x_i, b^y_i, b^z_i)^T
    b^+_i = b^x_i + im*b^y_i    ->     nb^+_i = nb^x_i + im*nb^y_i
procedure
    1.rewrite operator of ith site in global axis and perform symmetry, get nai
    2.rewrite operator of jth site in global axis, get aj
    3.get phase factor between nai and aj
"""
function get_phase_factor(sym, Ri, Rj, T; tor=1e-5)
    if T == 1
        nai = round.(sym.s*Ri*[1, 0, 0] + complex(0, 1)*sym.s*Ri*[0, 1, 0], digits=5)
    elseif T == -1
        nai = round.(sym.s*Ri*[1, 0, 0] - complex(0, 1)*sym.s*Ri*[0, 1, 0], digits=5)
    else
        error("T error")
    end
    aj = round.(Rj*[1, 0, 0] + complex(0, 1)*Rj*[0, 1, 0], digits=5)
    factor = [ ele for ele in nai ./ aj if !isequal(real(ele), NaN) && !isequal(imag(ele), NaN)]
    for i in 1:length(factor)-1
        if abs.(round.(factor[i]-factor[i+1], digits=5)) .< tor continue end
        println(nai)
        println(aj)
        println(factor)
        return nothing
    end
    return factor[1]
end

"""
get_rotation_matrix_from_local_axis_to_global_axis
"""
function get_rotation_matrix_from_local_axis_to_global_axis(cartesian_magmoms; tor=1e-5)::Vector{Matrix{Float64}}
    Rs = Vector{Matrix{Float64}}()
    for magmom in cartesian_magmoms
        r = sqrt(sum(magmom.^2))
        d = sqrt(sum(magmom[1:2].^2))
        theta = acos(magmom[3]/r)
        if d < tor
            phi = 0.0
        elseif magmom[2] > -tor
            phi = acos(magmom[1]/d)
        else
            phi = 2*pi-acos(magmom[1]/d)
        end
        R = Matrix{Float64}([
             cos(phi)*cos(theta) -sin(phi) cos(phi)*sin(theta)
             sin(phi)*cos(theta)  cos(phi) sin(phi)*sin(theta)
            -sin(theta)           0.0      cos(theta)
        ])
        push!(Rs, R)
    end
    return Rs
end

end#(module SystemGroupRepresentation)

################################################MagnonBandStructure###################################################
######################################################################################################################
module MagnonBandStructure
using ..SystemStructure
using ..SystemNeighbor
using ..SystemSymmetry

export get_atom_and_index_offset_of_site_with_finite_magmom

"""
get_atom_and_index_offset_of_site_with_finite_magmom(cell; tor=1e-3)
"""
function get_atom_and_index_offset_of_site_with_finite_magmom(cell; tor=1e-3)
    for i in 1:cell["size"]
        if any(abs.(round.(cell["magmom"][i], digits=3)) .> tor)
            return cell["atom"][i], i-1
        end
    end
    error()
end

end#(module MagnonBandStructure)


using .SystemStructure
using .SystemSite
using .SystemSymmetry
using .SystemNeighbor
using .kSpaceSymmetry
using .SystemGroupRepresentation
using .MagnonBandStructure
end #(module System)

#=
function minimize_cell(cell)
    cell_ori = deepcopy(cell)
    cell["cartesian_coord"] = [cell["lattice"]*coord for coord in cell["frac_coord"]]
    sites = get_full_sites(cell)
    eqsites_list = []
    #num = sum([1 for atom in cell["atom"] if atom == cell["atom"][0]])
    for i in 1:cell["size"]
        for eqsites in eqsites_list
            if i in eqsites @goto next_i end
        end
        push!(eqsites_list, [i])
        for j in 1:cell["size"]
            for eqsites in eqsites_list
                if j in eqsites @goto next_j end
            end
            if is_site_equivalent(i, j, cell)
                push!(eqsites_list[end], j)
            end
            @label next_j
        end
        @label next_i
    end
    eqsites_list
end
function get_primitive_cell(cell)
    trans_syms = get_translation_group_symmetries(cell)
    if trans_syms["size"] == 1

    else
       println(trans_syms["size"])
    end
end
=#
