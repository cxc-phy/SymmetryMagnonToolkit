module MathTool
using LinearAlgebra

export generate_list
export get_transformation_matrix
export composite_matrice
export get_direction_cosine_of_vector
export get_distance
export optimize_constraints
export is_vector_list_same
export is_vector_perpendicular
export is_vector_paralleled
export is_vector_antiparalleled
export get_direction_vector_of_line
export get_normal_vector_of_plane
export get_inverse_of_matrix
export get_all_permutations
export combine_permutations

"""
generate_list(idxs::Int64...) -> Vector{Vector{Int64}}
"""
function generate_list(rngs::Vector{UnitRange{Int64}})::Vector{Vector{Int64}}
    if length(rngs) == 1
        list = generate_list(rngs[1])
    else
        list1 = generate_list(rngs[1])
        list2 = generate_list(rngs[2:end])
        list = [[i; j] for i in list1 for j in list2]
    end
    return list
end
function generate_list(rngs::UnitRange{Int64})::Vector{Vector{Int64}}
    list = [[i] for i in rngs]
    return list
end
function generate_list(idxs::Vector{Int64})::Vector{Vector{Int64}}
    return generate_list(idxs...)
end
function generate_list(idxs::Int64...)::Vector{Vector{Int64}}
    if length(idxs) == 1
        list = [[i] for i in 0:idxs[1]]
    else
        list1 = generate_list(idxs[1])
        list2 = generate_list(idxs[2:end]...)
        list = [[i; j] for i in list1 for j in list2]
    end
    return list
end

"""
get_transformation_matrix(basis_list)
"""
function get_transformation_matrix(basis_list)
    dim = size(basis_list[1], 1)
    Tmat = zeros(Complex{Float64}, dim, dim)
    for i in 1:dim
        Tmat[:, i] = basis_list[i]
    end
    return Tmat
end

"""
composite_matrice(sym::T; order::Integer)::T where T<:Any
composite_matrice(syms::T...)::T where T<:Any
"""
function composite_matrice(sym::T; order::Integer)::T where T<:Any
    if order == 1
        return sym
    elseif order > 1
        return composite_matrice([sym for i = 1:order]...)
    else
        error("Error occurs at composite_symmetries: order should be a positive integer")
    end
end
function composite_matrice(syms::T...)::T where T<:Any
    if length(syms) == 1
        return syms[1]
    elseif length(syms) == 2
        return syms[1]*syms[2]
    else
        return composite_matrice(composite_matrice(syms[1], syms[2]), syms[3:end]...)
    end
end

"""
get_direction_cosine_of_vector(vec::Vector{Float64}) -> Vector{Float64}
"""
function get_direction_cosine_of_vector(vec::Vector{Float64}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64})::Vector{Float64}
    x = x/norm(x)
    y = y/norm(y)
    z = z/norm(z)
    vec = vec/norm(vec)
    dir_cosine = round.([sum(x.*vec), sum(y.*vec), sum(z.*vec)], digits=3)
    return dir_cosine
end

"""
get_distance(coord1, coord2) -> Float64
"""
function get_distance(coord1, coord2)::Float64
    dist = round.(norm(coord1-coord2), digits=5)
    return dist
end

"""
optimize_constraints(constrs; tor=1e-3)
"""
function optimize_constraints(constrs; tor=1e-3)
    #=
    constrs:store constraints of k
    size(constrs)=(3,4) ->
        constrs[1,1]*u+constrs[1,2]*v+constrs[1,3]*w+constrs[1,4]=0
        constrs[2,1]*u+constrs[2,2]*v+constrs[2,3]*w+constrs[1,4]=0
        constrs[3,1]*u+constrs[3,2]*v+constrs[3,3]*w+constrs[1,4]=0
    nconstrs_list:store several reduced constraints of k
    nconstrs:store one reduced constraints of k
    size(nconstrs)=(3,4) ->
        constrs[1,1]*u+constrs[1,4]=0
        constrs[2,1]*u+constrs[2,2]*v+constrs[1,4]=0
        constrs[3,1]*u+constrs[3,2]*v+constrs[3,3]*w+constrs[1,4]=0
    used_row: store the rowth row that has already been reduced
    =#
    for row in [1, 2, 3]
        if all(abs.(round.(constrs[row, 1:3], digits=3)) .< tor) && abs(round(constrs[row, 4], digits=3)) > tor
            @goto no_result
        end
    end
    nconstrs_list = []
    used_row = []
    for var in [3, 2, 1]
        for row in [3, 2, 1]
            if abs(round(constrs[row, var], digits=3)) < tor || row in used_row
                #=
                condition: the rowth row gives no constraints on varth variable or is already used
                statement: continue
                =#
                continue
            end
            if all([i==var ? 0 : abs(round(constrs[row, i], digits=3)) for i in 1:3] .< tor)
                #=
                condition: the rowth row in constrs can give exact values of the varth variable
                statement:
                #1.put the rowth row into used_row
                #2.get exact values
                #3.go through other rows not in used_row, check if each value is compatible with them
                    if true,
                        (1)store value in kvarsols
                        (2)if the row gives the same result as current one, put it into used_row
                    if false, do nothing
                #4.append rows in used_idx into used_row
                kvarsols: store solution for the varth variable
                used_idx: store the idxth row giving the same constraints as the rowth row
                =#
                push!(used_row, row)
                used_idx = []
                sol = -constrs[row, 4]/constrs[row, var]
                kvarsol = [[i == var ? -1 : 0 for i = 1:3];sol]
                for idx in [3, 2, 1]
                    if abs(round(constrs[idx, var], digits=3)) < tor || idx in used_row continue end
                    constr = constrs[idx, :] - constrs[idx, var]/kvarsol[var]*kvarsol[:]
                    if any(abs.(round.(constr[1:3], digits=3)) .> tor)
                        #constrs[idx, :] .-= constrs[idx, var]/kvarsol[var]*kvarsol[:]
                        continue
                    end
                    if !(idx in used_idx) push!(used_idx, idx) end
                    if abs(mod(round(constr[4], digits=3), 1)) < tor continue end
                    #if abs(round(constr[4], digits=3)) < tor continue end
                    @goto no_result
                end
                kvarsols=[kvarsol]
                append!(used_row, used_idx)
            else
                #=
                condition: the rowth row in constrs gives constraint on the varth variable
                statement:
                #1.check if the coefficent of the varth variable (constrs[row, var]) is +1/-1
                    if true, store constraint in kvarsols
                    if false and exist another row having true result, do nothing
                    if false and exist no other row having true result, store constraint in kvarsols
                #2.if constraint is stored, eliminate the varth variable in the remaining constraints
                *the purpose is to avoid division due to the constraints is actually mod 1
                kvarsols: store solution for the varth variable
                =#
                if abs(round(abs(constrs[row, var])-1, digits=3)) > tor
                    @goto next_row
                    for idx = 1:row-1
                        if abs(round(abs(constrs[idx, var])-1, digits=3)) < tor @goto next_row end
                    end
                end
                push!(used_row, row)
                kvarsols = []
                push!(kvarsols, constrs[row, :]*-constrs[row, var]/abs(constrs[row, var]))
                for idx in [3, 2, 1]
                    if abs(round(constrs[idx, var], digits=3)) < tor || idx in used_row continue end
                    constrs[idx, :] .-= constrs[idx, var]/constrs[row, var]*constrs[row, :]
                end
            end
            #=

            add kvarsols to nconstrs_list
            =#
            if nconstrs_list == []
                nconstrs_list = [kvarsol for kvarsol in kvarsols]
            else
                nconstrs_list = [[kvarsol; nconstrs] for kvarsol in kvarsols for nconstrs in nconstrs_list]
            end
            break
            @label next_row
        end
        #=the varth variable has no exact values or constraints=#
        if nconstrs_list == []
            nconstrs_list = [Float64.([0, 0, 0, 0])]
        elseif length(nconstrs_list[1]) == 4*(3-var)
            nconstrs_list = [[Float64.([0, 0, 0, 0]); nconstrs] for nconstrs in nconstrs_list]
        elseif length(nconstrs_list[1]) != 4*(4-var)
            error("Error occurs at get_constraints_list_of_symmetry")
        end
    end
    nconstrs_list = [[ele[4*(i-1)+j] for i=1:3, j=1:4] for ele in nconstrs_list]
    for var in [1, 2]
        for i in eachindex(nconstrs_list)
            nconstrs = nconstrs_list[i]
            if abs(round(nconstrs[var, var], digits=3)) < tor continue end
            nconstrs = [ idx in var+1:3 ? nconstrs[idx, :] .-= nconstrs[idx, var]/nconstrs[var, var]*nconstrs[var, :] : nconstrs[idx, :] for idx = 1:3]
            nconstrs_list[i] = [nconstrs[i][j] for i=1:3, j=1:4]
        end
    end
    return nconstrs_list

    @label no_result
    return nothing
end

function is_vector_list_same(vecs1, vecs2; tor=1e-3)::Bool
    for vec1 in vecs1
        for vec2 in vecs2
            if all(abs.(round.(vec1-vec2, digits=3)) .< tor) @goto next_vec end
        end
        return false
        @label next_vec
    end
    return true
end

function is_vector_perpendicular(vec1, vec2; tor=1e-3)::Bool
    if all(abs.(round.(dot(vec1, vec2), digits=3)) .< tor)
        return true
    else
        return false
    end
end

function is_vector_paralleled(vec1, vec2; tor=1e-3)::Bool
    if all(abs.(round.(cross(vec1, vec2), digits=3)) .< tor) && dot(vec1, vec2) > tor
        return true
    else
        return false
    end
end

function is_vector_antiparalleled(vec1, vec2; tor=1e-3)::Bool
    if all(abs.(round.(cross(vec1, vec2), digits=3)) .< tor) && dot(vec1, vec2) < tor
        return true
    else
        return false
    end
end

function get_direction_vector_of_line(line::Vector{Float64}; tor=1e-3)
    nline = zeros(Float64, 3, 4)
    for i in 1:3
        if abs(round(line[i], digits=3)) < tor continue end
        if i == 1
            nline[1, 1] = 1.0
        else
            if abs(round(line[i-1], digits=3)) < tor
                nline[i, 1] = 1.0
            else
                nline[i, 1] = line[i]/line[i-1]
            end
        end
    end
    return get_direction_vector_of_line(nline)
end
function get_direction_vector_of_line(line::Matrix{Float64}; tor=1e-3)
    uvwlist = [
        [0.0, 0.0, 0.0, 1],
        [0.0, 0.0, 0.5, 1],
        [0.0, 0.5, 0.0, 1],
        [0.0, 0.5, 0.5, 1],
        [0.5, 0.0, 0.0, 1],
        [0.5, 0.0, 0.5, 1],
        [0.5, 0.5, 0.0, 1],
        [0.5, 0.5, 0.5, 1]
    ]
    ptstart = line[:, 4]
    ptend = nothing
    for uvw in uvwlist
        ptend = line*uvw
        if round(norm(ptstart-ptend),digits=3) > tor break end
    end
    if ptend == nothing
        return nothing
    end
    dirvec = (ptend-ptstart)/norm(ptend-ptstart)
    return dirvec
end

"""
#size(kp) = (3, 4) -> (
    kx = kp[1,1]*u+kp[1,2]*v+kp[1,3]*w+kp[1,4],
    ky = kp[2,1]*u+kp[2,2]*v+kp[2,3]*w+kp[2,4],
    kz = kp[3,1]*u+kp[3,2]*v+kp[3,3]*w+kp[3,4]
)
"""
function get_normal_vector_of_plane(plane; tor=1e-3)
    uvwlist = [
        [0.0, 0.0, 0.0, 1],
        [0.0, 0.0, 0.5, 1],
        [0.0, 0.5, 0.0, 1],
        [0.0, 0.5, 0.5, 1],
        [0.5, 0.0, 0.0, 1],
        [0.5, 0.0, 0.5, 1],
        [0.5, 0.5, 0.0, 1],
        [0.5, 0.5, 0.5, 1]
    ]
    ptstart = plane[:, 4]
    #=find one line on plane=#
    ptend1 = nothing
    for uvw in uvwlist
        ptend1 = plane*uvw
        if round(norm(ptstart-ptend1), digits=3) > tor break end
    end
    if ptend1 == nothing
        return nothing
    end
    dirvec1 = (ptend1-ptstart)/norm(ptend1-ptstart)
    #=find another line on each plane=#
    ptend2 = nothing
    for uvw in uvwlist
        ptend2 = plane*uvw
        if round(norm(ptstart-ptend2),digits=3) > tor &&
           round(norm(ptend1-ptend2),digits=3) > tor
            break
        end
    end
    if ptend2 == nothing
        return nothing
    end
    dirvec2 = (ptend2-ptstart)/norm(ptend2-ptstart)
    #=get normal vectors=#
    normvec = cross(dirvec1, dirvec2)
    return normvec
end

"""
get_inverse_of_matrix(sym::Matrix{Int64}; tor=1e-3) -> Matrix{Int64}
"""
function get_inverse_of_matrix(sym::Matrix{T}; tor=1e-3)::Matrix{T} where {T<:Union{Int64, Float64}}
    sym_inv_prime = inv(sym)
    if T == Int64
        sym_inv = Int64.(round.(sym_inv_prime, digits=1))
    elseif T == Float64
        sym_inv = round.(sym_inv_prime, digits=1)
    end
    if all(sym_inv-sym_inv_prime .< tor)
        return sym_inv
    else
        error("Error occurs at get_inverse_of_matrix:inverse of matrix can not be found")
    end
end

"""
get_all_permutations(lists::Vector{Vector{Int64}}) -> Vector{Vector{Int64}}
get_all_permutations(lists::Vector{Int64}) -> Vector{Vector{Int64}}
"""
function get_all_permutations(lists::Vector{Vector{Int64}})
    perm_previous = []
    for list in lists
        perm_latter = get_all_permutations(list)
        if perm_previous == []
            perm_previous = perm_latter
        else
            perm_previous = [[perm_previous[i];perm_latter[j]] for i=1:length(perm_previous) for j=1:length(perm_latter)]
        end
    end
    permutations = perm_previous
    return permutations
end

#=
function get_all_permutations(list::Vector{Int64})
    if length(list) > 8
        error("Error occurs at get_all_permutations:can not get permutations for order >10")
    elseif length(list) == 1
        return [[list[1]]]
    else
        permutations = []
        for i in eachindex(list)
            x = list[i]
            if i == 1
                nlist = list[2:end]
            elseif i == length(list)
                nlist = list[1:end-1]
            else
                nlist = [list[1:i-1];list[i+1:end]]
            end
            for ele in get_all_permutations(nlist)
                push!(permutations,[x;ele])
            end
        end
        return permutations
    end
end
=#
function get_all_permutations(list::Vector{Int64})
    function backtrack(first=1)
        if first == num
            push!(result, list[:])
        end
        for i in first:num
            list[first], list[i] = list[i], list[first]
            backtrack(first+1)
            list[first], list[i] = list[i], list[first]
        end
    end
    num=length(list)
    result=Vector{Vector{Int64}}()
    backtrack()
    return result
end


function combine_permutations(perms...)
    if length(perms) == 1
        perm = perms[1]
        return perm
    elseif length(perms) == 2
        perm = [perms[2][i] for i in perms[1]]
        return perm
    else
        return combine_permutations(combine_permutations(perms[1], perms[2]), perms[3:end]...)
    end
end

end
