try
    typeof(kSpace)
catch UndefVarError
    include("../representation-tool/kspace.jl")
end
try
    typeof(PointSymmetry)
catch UndefVarError
    include("../symmetry-operation-tool/pointsymmetry.jl")
end
try
    typeof(GroupWithMultiplicationTable)
catch UndefVarError
    include("../representation-tool/groupwithmultiplicationtable.jl")
end
######################################################################################################################
######################################################PointGroup######################################################
######################################################################################################################
module PointGroup
include("../database/pg_database.jl")

####################################################PointSymmetry#####################################################
######################################################################################################################
module PointSymmetry
using LinearAlgebra
using ...PointSymmetry
using ...MathTool
using ..PointGroupDataBase

export get_symmetry_order
export get_isometry_type
export get_isometry_axis
export get_proper_part
export get_sense_of_rotation
export is_sym1_commuting_with_sym2
export is_sym_commuting_with_syms
export composite_symmetries
export conjugate_symmetry

"""
get_symmetry_order(encoded::pg_encodetype) -> Int64
get_symmetry_order(sym::pg_symtype) -> Int64
#get order n of a symmetry operation, sym^n=I
"""
function get_symmetry_order(encoded::pg_encodetype)::Int64
    sym = pg_decode_symmetry(encoded)
    return get_symmetry_order_of_point_symmetry(sym)
end
function get_symmetry_order(sym::pg_symtype)::Int64
    return get_symmetry_order_of_point_symmetry(sym)
end

"""
get_isometry_type(encoded::pg_encodetype) -> Int64
get_isometry_type(sym::pg_symtype) -> Int64
#get isometry type of a symmetry operation
    operation    -6 -4 -3 -2 -1  1  2  3  4  6
    index         1  2  3  4  5  6  7  8  9 10
"""
function get_isometry_type(encoded::pg_encodetype)::Int64
    sym = pg_decode_symmetry(encoded)
    return get_isometry_type_of_point_symmetry(sym)
end
function get_isometry_type(sym::pg_symtype)::Int64
    return get_isometry_type_of_point_symmetry(sym)
end

"""
get_isometry_axis(encoded::pg_encodetype) -> Vector{Int64}
get_isometry_axis(sym::pg_symtype) -> Vector{Int64}
#get isometry axis of a symmetry operation
"""
function get_isometry_axis(encoded::pg_encodetype)::Vector{Int64}
    sym = pg_decode_symmetry(encoded)
    return get_isometry_axis_of_point_symmetry(sym)
end
function get_isometry_axis(sym::pg_symtype)::Vector{Int64}
    return get_isometry_axis_of_point_symmetry(sym)
end

"""
get_proper_part(encoded::pg_encodetype; encode::Bool) -> Union{pg_encodetype, pg_symtype}
get_proper_part(sym::pg_symtype; encode::Bool) -> Union{pg_encodetype, pg_symtype}
#get proper part of a symmetry operation
"""
function get_proper_part(encoded::pg_encodetype; encode::Bool)::Union{pg_encodetype, pg_symtype}
    sym = pg_decode_symmetry(encoded)
    return get_proper_part(sym, encode=encode)
end
function get_proper_part(sym::pg_symtype; encode::Bool)::Union{pg_encodetype, pg_symtype}
    proper = get_proper_part_of_point_symmetry(sym)
    if encode
        return pg_encode_symmetry(proper)
    else
        return proper
    end
end

"""
get_sense_of_rotation(encoded::pg_encodetype) -> Int64
get_sense_of_rotation(sym::pg_symtype) -> Int64
#get sense of rotation of a symmetry operation
"""
function get_sense_of_rotation(encoded::pg_encodetype)::Int64
    sym = pg_decode_symmetry(encoded)
    return get_sense_of_rotation_of_point_symmetry(sym)
end
function get_sense_of_rotation(sym::pg_symtype)::Int64
    return get_sense_of_rotation_of_point_symmetry(sym)
end

"""
is_sym1_commuting_with_sym2(;sym1::T, sym2::T) -> Bool where T<:Union{pg_encodetype, pg_symtype}
#check if sym1 and sym2 commute with each other
"""
function is_sym1_commuting_with_sym2(;sym1::T, sym2::T)::Bool where T<:Union{pg_encodetype, pg_symtype}
    sym12 = composite_symmetries(sym1, sym2, encode=false)
    sym21 = composite_symmetries(sym2, sym1, encode=false)
    if all(sym12 .== sym21)
        return true
    else
        return false
    end
end
"""
is_sym_commuting_with_syms(encoded::pg_encodetype, encodeds::Vector{pg_encodetype}) -> Bool
is_sym_commuting_with_syms(sym::pg_symtype, syms::Vector{pg_symtype}) -> Bool
#check if sym commutes with all elements in syms
"""
function is_sym_commuting_with_syms(encoded::pg_encodetype, encodeds::Vector{pg_encodetype})::Bool
    sym = pg_decode_symmetry(encoded)
    syms = [pg_decode_symmetry(encoded) for encoded in encodeds]
    return is_sym_commuting_with_syms(sym, syms)
end
function is_sym_commuting_with_syms(sym::pg_symtype, syms::Vector{pg_symtype})::Bool
    for ele in syms
        if !is_sym1_commuting_with_sym2(sym1=sym, sym2=ele) return false end
    end
    return true
end

"""
composite_symmetries(encoded::pg_encodetype; order::Integer, encode::Bool) -> Union{pg_encodetype, pg_symtype}
composite_symmetries(sym::pg_symtype; order::Integer, encode::Bool)::Union{pg_encodetype, pg_symtype}
composite_symmetries(encodeds::pg_encodetype...; encode::Bool)::Union{pg_encodetype, pg_symtype}
composite_symmetries(syms::pg_symtype...; encode::Bool)::Union{pg_encodetype, pg_symtype}
"""
function composite_symmetries(encoded::pg_encodetype; order::Integer, encode::Bool)::Union{pg_encodetype, pg_symtype}
    sym = pg_decode_symmetry(encoded)
    return composite_symmetries(sym, order=order, encode=encode)
end
function composite_symmetries(sym::pg_symtype; order::Integer, encode::Bool)::Union{pg_encodetype, pg_symtype}
    if order == 1
        if encode
            return pg_encode_symmetry(sym)
        else
            return sym
        end
    elseif order > 1
        return composite_symmetries([sym for i = 1:order]..., encode=encode)
    else
        error("Error occurs at composite_symmetries: order should be a positive integer")
    end
end
function composite_symmetries(encodeds::pg_encodetype...; encode::Bool)::Union{pg_encodetype, pg_symtype}
    syms = (pg_decode_symmetry(encoded) for encoded in encodeds)
    return composite_symmetries(syms..., encode=encode)
end
function composite_symmetries(syms::pg_symtype...; encode::Bool)::Union{pg_encodetype, pg_symtype}
    if length(syms) == 1
        if encode
            return pg_encode_symmetry(syms[1])
        else
            return syms[1]
        end
    elseif length(syms) == 2
        if encode
            return pg_encode_symmetry(syms[1]*syms[2])
        else
            return syms[1]*syms[2]
        end
    else
        return composite_symmetries(composite_symmetries(syms[1], syms[2], encode=false), syms[3:end]..., encode=encode)
    end
end

"""
conjugate_symmetry(sym_encoded::pg_encodetype, conjg_sym_encoded::pg_encodetype; encode::Bool) -> Union{pg_encodetype, pg_symtype}
conjugate_symmetry(sym::pg_symtype, conjg_sym::pg_symtype; encode::Bool) -> Union{pg_encodetype, pg_symtype}
#matrices are in the same coordinate systems by default
#conjg_sym^{-1}*sym*conjg_sym
"""
function conjugate_symmetry(sym_encoded::pg_encodetype, conjg_sym_encoded::pg_encodetype; encode::Bool)::Union{pg_encodetype, pg_symtype}
    sym = pg_decode_symmetry(sym_encoded)
    conjg_sym = pg_decode_symmetry(conjg_sym_encoded)
    return conjugate_symmetry(sym, conjg_sym, encode=encode)
end
function conjugate_symmetry(sym::pg_symtype, conjg_sym::pg_symtype; encode::Bool)::Union{pg_encodetype, pg_symtype}
    return composite_symmetries(Int64.(inv(conjg_sym)), sym, conjg_sym, encode=encode)
end

end#(module PointSymmetry)

###################################################GroupInformation###################################################
######################################################################################################################
module GroupInformation
using LinearAlgebra
using ...MathTool
using ..PointGroupDataBase
using ..PointSymmetry
import ...GroupWithMultiplicationTable as gwmt

export get_isometry_type_of_conjugacy_class
export get_full_conjugacy_classes
export get_holohedry_of_group
export isa_subgroup
export is_group_same
export get_permutation_to_standard_point_group
export get_isomorphism_to_standard_point_group
export get_multiplication_table
export conjugate_group
export get_class_of_point_group
export get_notation_of_point_group
export get_isometry_type_table
export is_group_abelian
export isa_set
export isa_group
export get_full_symmetry_operations
export generate_full_symmetry_operations
export get_full_subgroups
export get_certain_subgroup

#=
get_isometry_type_of_conjugacy_class(conjg_class::Vector{pg_encodetype}) -> Int64
get_isometry_type_of_conjugacy_class(conjg_class::Vector{pg_symtype}) -> Int64
get_full_conjugacy_classes(pgnum::Int64; encode::Bool) -> Union{Vector{Vector{pg_encodetype}}, Vector{Vector{pg_symtype}}}
get_full_conjugacy_classes(syms::Vector{pg_symtype}; encode::Bool) -> Union{Vector{Vector{pg_encodetype}}, Vector{Vector{pg_symtype}}}
get_full_conjugacy_classes(syms::Vector{pg_encodetype}; encode::Bool) -> Union{Vector{Vector{pg_encodetype}}, Vector{Vector{pg_symtype}}}
=#
#=
isa_subgroup(group1::Vector{pg_encodetype}, group2::Vector{pg_encodetype}) -> Bool
is_group_same(group1::Vector{Vector{pg_symtype}}, group2::Vector{Vector{pg_encodetype}}) -> Bool
is_group_same(group1::Vector{Vector{pg_encodetype}}, group2::Vector{Vector{pg_symtype}}) -> Bool
is_group_same(group1::Vector{Vector{pg_encodetype}}, group2::Vector{Vector{pg_encodetype}}) -> Bool
get_permutation_to_standard_point_group(group::Vector{pg_symtype}; std_group::Vector{pg_symtype}) -> Vector{Int64}
get_permutation_to_standard_point_group(group::Vector{pg_symtype}; std_group::Union{Vector{pg_encodetype}, Nothing}=nothing) -> Vector{Int64}
get_permutation_to_standard_point_group(group::Vector{pg_encodetype}; std_group::Union{Vector{pg_encodetype}, Nothing}=nothing) -> Vector{Int64}
get_isomorphism_to_standard_point_group(group::Vector{pg_symtype}; std_group::Vector{pg_symtype}) -> Vector{Int64}
get_isomorphism_to_standard_point_group(group::Vector{pg_symtype}; std_group::Union{Vector{pg_encodetype}, Nothing}=nothing) -> Vector{Int64}
get_isomorphism_to_standard_point_group(group::Vector{pg_encodetype}; std_group::Union{Vector{pg_encodetype}, Nothing}=nothing) -> Vector{Int64}
permute_multiplication_table(permutation, multable)
get_multiplication_table(pgnum::Int64) -> Matrix{Int64}
get_multiplication_table(encodeds::Vector{pg_encodetype}) -> Matrix{Int64}
get_multiplication_table(syms::Vector{pg_symtype}) -> Matrix{Int64}
conjugate_group(pgnum::Int64, conjg_sym_encoded::pg_encodetype; encode::Bool) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
conjugate_group(group::Vector{pg_symtype}, conjg_sym::pg_symtype; encode::Bool) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
get_class_of_point_group(syms::Union{pg_encodetype, pg_symtype}...) -> Int64
get_notation_of_point_group(syms::Union{Vector{pg_encodetype}, Vector{pg_symtype}})
get_notation_of_point_group(idx::Int64)
get_isometry_type_table(syms::Union{Vector{pg_encodetype}, Vector{pg_symtype}}) -> Vector{Int64}
get_isometry_type_table(encodeds::pg_encodetype...) -> Vector{Int64}
get_isometry_type_table(syms::pg_symtype...) -> Vector{Int64}
is_group_abelian(pgnum::Int64) -> Bool
isa_set(list::Union{Vector{pg_encodetype}, Vector{pg_symtype}}) -> Bool
isa_set(encodeds::pg_encodetype...) -> Bool
isa_set(list::pg_symtype...) -> Bool
isa_group(list::Union{Vector{pg_encodetype}, Vector{pg_symtype}}) -> Bool
isa_group(list::pg_encodetype...) -> Bool
isa_group(list::pg_symtype...) -> Bool
get_full_symmetry_operations(pgnum::Int64; encode::Bool) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
generate_full_symmetry_operations(syms::Union{Vector{pg_encodetype}, Vector{pg_symtype}}; encode::Bool, sort::Bool=false) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
generate_full_symmetry_operations(encodeds::pg_encodetype...; encode::Bool, sort::Bool=false) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
generate_full_symmetry_operations(syms::pg_symtype...; encode::Bool, sort::Bool=false) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
get_full_subgroups(pgnum::Int64; encode::Bool) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
get_certain_subgroup(pgnum::Int64, subindx::Int64; encode::Bool) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
=#

"""
get_isometry_type_of_conjugacy_class(conjg_class::Vector{pg_encodetype}) -> Int64
get_isometry_type_of_conjugacy_class(conjg_class::Vector{pg_symtype}) -> Int64
isometry_type_table
index      1  2  3  4  5  6  7  8  9 10
operation -6 -4 -3 -2 -1  1  2  3  4  6
"""
function get_isometry_type_of_conjugacy_class(conjg_class::Vector{pg_encodetype})::Int64
    conjg_class = [pg_decode_symmetry(sym) for sym in conjg_class]
    return get_isometry_type_of_conjugacy_class(conjg_class)
end
function get_isometry_type_of_conjugacy_class(conjg_class::Vector{pg_symtype})::Int64
    iso_type = 0
    for sym in conjg_class
        if iso_type == 0
            iso_type = get_isometry_type(sym)
        elseif iso_type != get_isometry_type(sym)
            error("Error occurs at get_isometry_type_of_conjugacy_class: symmetries in the same conjugacy class should have same isometry type")
        end
    end
    if iso_type == 0 error("Error occurs at get_isometry_type_of_conjugacy_class: can not find isometry type of conjugacy class") end
    return iso_type
end

"""
get_full_conjugacy_classes(pgnum::Int64; encode::Bool) -> Union{Vector{Vector{pg_encodetype}}, Vector{Vector{pg_symtype}}}
get_full_conjugacy_classes(syms::Vector{pg_symtype}; encode::Bool) -> Union{Vector{Vector{pg_encodetype}}, Vector{Vector{pg_symtype}}}
get_full_conjugacy_classes(syms::Vector{pg_encodetype}; encode::Bool) -> Union{Vector{Vector{pg_encodetype}}, Vector{Vector{pg_symtype}}}
sort_conjg_classes(conjg_classes::Vector{Vector{pg_encodetype}}; encode::Bool) -> Union{Vector{Vector{pg_encodetype}}, Vector{Vector{pg_symtype}}}
#sort conjugacy classes in arrangement of operation[1, -6, -4, -3, -2, -1, 2, 3, 4, 6] or index[6, 1, 2, 3, 4, 5, 7, 8, 9]
#isometry_type_table
 index      1  2  3  4  5  6  7  8  9 10
 operation -6 -4 -3 -2 -1  1  2  3  4  6
"""
function get_full_conjugacy_classes(pgnum::Int64; encode::Bool)::Union{Vector{Vector{pg_encodetype}}, Vector{Vector{pg_symtype}}}
    syms = get_full_symmetry_operations(pgnum, encode=true)
    return get_full_conjugacy_classes(syms, encode=encode)
end
function get_full_conjugacy_classes(syms::Vector{pg_symtype}; encode::Bool)::Union{Vector{Vector{pg_encodetype}}, Vector{Vector{pg_symtype}}}
    encodeds = [pg_encode_symmetry(sym) for sym in syms]
    return get_full_conjugacy_classes(encodeds, encode=encode)
end
function get_full_conjugacy_classes(syms::Vector{pg_encodetype}; encode::Bool)::Union{Vector{Vector{pg_encodetype}}, Vector{Vector{pg_symtype}}}
    syms_remain = [sym for sym in syms if sym !== 16484]
    conjg_classes = Vector{Vector{Int64}}([[16484]])
    while syms_remain != []
        class = Vector{Int64}()
        push!(class, syms_remain[1])
        for sym in syms
            nsym = conjugate_symmetry(syms_remain[1], sym, encode=true)
            if !(nsym in syms) error("Error occurs at get_conjugacy_classes_of_point_groups: group not closed") end
            if !(nsym in class) push!(class, nsym) end
        end
        push!(conjg_classes, class)
        syms_remain = [sym for sym in syms_remain if !(sym in class)]
    end
    return sort_conjg_classes(conjg_classes; encode=encode)
end
function sort_conjg_classes(conjg_classes::Vector{Vector{pg_encodetype}}; encode::Bool)::Union{Vector{Vector{pg_encodetype}}, Vector{Vector{pg_symtype}}}
    iso_types = []
    for class in conjg_classes
        push!(iso_types, get_isometry_type_of_conjugacy_class(class))
    end
    iso_types = [6; deleteat!(sort(iso_types), findall(x->x==6, sort(iso_types)))]
    sorted_conjg_classes = Vector{Vector{Int64}}()
    for iso_type in iso_types
        for i in eachindex(conjg_classes)
            class = conjg_classes[i]
            if get_isometry_type_of_conjugacy_class(class)==iso_type
                push!(sorted_conjg_classes, class)
                deleteat!(conjg_classes, i)
                break
            end
        end
    end
    if !encode sorted_conjg_classes = [[pg_decode_symmetry(encoded) for encoded in class] for class in sorted_conjg_classes] end
    return sorted_conjg_classes
end

"""
get_holohedry_of_group(idx)
"""
function get_holohedry_of_group(idx)
    return pg_data[idx].holohedry
end

"""
isa_subgroup(subnum::Int64, groupnum::Int64) -> Bool
isa_subgroup(group1::Vector{pg_encodetype}, group2::Vector{pg_encodetype}) -> Bool
"""
function isa_subgroup(subnum::Int64, groupnum::Int64)::Bool
    group1 = get_full_symmetry_operations(subnum, encode=true)
    group2 = get_full_symmetry_operations(groupnum, encode=true)
    return isa_subgroup(group1, group2)
end
function isa_subgroup(group1::Vector{pg_encodetype}, group2::Vector{pg_encodetype})::Bool
    for encoded in group1
        if !(encoded in group2) return false end
    end
    return true
end


"""
is_group_same(group1::Vector{pg_symtype}, group2::Vector{pg_encodetype}) -> Bool
is_group_same(group1::Vector{pg_encodetype}, group2::Vector{pg_symtype}) -> Bool
is_group_same(group1::Vector{pg_encodetype}, group2::Vector{pg_encodetype}) -> Bool
#check if two groups are exactly the same without any isomorphism
"""
function is_group_same(group1::Vector{pg_symtype}, group2::Vector{pg_encodetype})::Bool
    group1 = [pg_encode_symmetry(sym) for sym in group1]
    return is_group_same(group1, group2)
end
function is_group_same(group1::Vector{pg_encodetype}, group2::Vector{pg_symtype})::Bool
    group2 = [pg_encode_symmetry(sym) for sym in group2]
    return is_group_same(group1, group2)
end
function is_group_same(group1::Vector{pg_encodetype}, group2::Vector{pg_encodetype})::Bool
    if length(group1) != length(group2) return false end
    for sym in group1
        if !(sym in group2) return false end
    end
    return true
end

"""
get_permutation_to_standard_point_group(group::Vector{pg_symtype}; std_group::Vector{pg_symtype}) -> Vector{Int64}
get_permutation_to_standard_point_group(group::Vector{pg_symtype}; std_group::Union{Vector{pg_encodetype}, Nothing}=nothing) -> Vector{Int64}
get_permutation_to_standard_point_group(group::Vector{pg_encodetype}; std_group::Union{Vector{pg_encodetype}, Nothing}=nothing) -> Vector{Int64}
#find correspondence between group and standard point group when they are exactly the same
"""
function get_permutation_to_standard_point_group(group::Vector{pg_symtype}; std_group::Vector{pg_symtype})::Vector{Int64}
    group = [pg_encode_symmetry(sym) for sym in group]
    std_group = [pg_encode_symmetry(sym) for sym in group]
    return get_permutation_to_standard_point_group(group, std_group=std_group)
end
function get_permutation_to_standard_point_group(group::Vector{pg_symtype}; std_group::Union{Vector{pg_encodetype}, Nothing}=nothing)::Vector{Int64}
    group = [pg_encode_symmetry(sym) for sym in group]
    return get_permutation_to_standard_point_group(group, std_group=std_group)
end
function get_permutation_to_standard_point_group(group::Vector{pg_encodetype}; std_group::Union{Vector{pg_encodetype}, Nothing}=nothing)::Vector{Int64}
    if std_group == nothing
        std_pgnum = get_class_of_point_group(group...)
        std_group = get_full_symmetry_operations(std_pgnum, encode=true)
    end
    permutation = []
    for encoded in group
        for i in eachindex(std_group)
            if encoded == std_group[i]
                append!(permutation, i)
                break
            end
        end
    end
    if length(permutation) != length(group)
        error("Error occurs at get_permutation_to_standard_point_group: group and standard group are not same, please use get_isomorphism_to_standard_point_group instead")
    end
    return permutation
end

"""
get_isomorphism_to_standard_point_group(group::Vector{pg_symtype}; std_group::Vector{pg_symtype}) -> Vector{Int64}
get_isomorphism_to_standard_point_group(group::Vector{pg_symtype}; std_group::Union{Vector{pg_encodetype}, Nothing}=nothing) -> Vector{Int64}
get_isomorphism_to_standard_point_group(group::Vector{pg_encodetype}; std_group::Union{Vector{pg_encodetype}, Nothing}=nothing) -> Vector{Int64}
#find correspondence between group and standard point group when they are isomorphic but not exactly the same
#symmetry operations of standard point group are in arrangement of isometry types = [6, 1, 2, 3, 4, 5, 7, 8, 9, 10]
"""
function get_isomorphism_to_standard_point_group(group::Vector{pg_symtype}; std_group::Vector{pg_symtype})::Vector{Int64}
    group = [pg_encode_symmetry(sym) for sym in group]
    std_group = [pg_encode_symmetry(sym) for sym in group]
    return get_isomorphism_to_standard_point_group(group, std_group=std_group)
end
function get_isomorphism_to_standard_point_group(group::Vector{pg_symtype}; std_group::Union{Vector{pg_encodetype}, Nothing}=nothing)::Vector{Int64}
    group = [pg_encode_symmetry(sym) for sym in group]
    return get_isomorphism_to_standard_point_group(group, std_group=std_group)
end
function get_isomorphism_to_standard_point_group(group::Vector{pg_encodetype}; std_group::Union{Vector{pg_encodetype}, Nothing}=nothing)::Vector{Int64}
    if std_group == nothing
        #get symmetry opertaions of standard point group
        std_pgnum = get_class_of_point_group(group...)
        std_group = get_full_symmetry_operations(std_pgnum, encode=true)
    end
    if length(group) != length(std_group) error("error") end

    #####sort group in arrangement of isometry types = [6, 1, 2, 3, 4, 5, 7, 8, 9, 10]#####
    sorted_conjg_classes = get_full_conjugacy_classes(group, encode=true)
    sorted_group = Vector{Int64}()
    for class in sorted_conjg_classes
        append!(sorted_group, class)
    end

    #####get isomorphism from group to sorted_group#####
    perm1 = get_permutation_to_standard_point_group(group, std_group=sorted_group)

    #####get isomorphism from sorted_group to std_group#####
    multable = get_multiplication_table(sorted_group)
    std_multable = get_multiplication_table(std_group)
    lists = Vector{Vector{Int64}}()
    iso_type_class = []
    iso_type = 0
    for class in sorted_conjg_classes
        new_iso_type = get_isometry_type_of_conjugacy_class(class)
        if new_iso_type != iso_type
            push!(iso_type_class, class)
            iso_type = new_iso_type
        else
            append!(iso_type_class[end], class)
        end
    end
    for i = 1:length(iso_type_class)
        if i == 1
            list = [j for j=1:length(iso_type_class[i])]
        else
            start = sum(length(iso_type_class[k]) for k = 1:i-1)
            finish = start+length(iso_type_class[i])
            list = [j for j=start+1:finish]
        end
        push!(lists,list)
    end
    if *([factorial(length(list)) for list in lists]...) > 1e7
        println("number of test permutations is too large")
        #=
        for i = 1:length(sorted_conjg_classes)
            if i == 1
                list = [j for j=1:length(sorted_conjg_classes[i])]
            else
                start = sum(length(sorted_conjg_classes[k]) for k = 1:i-1)
                finish = start+length(sorted_conjg_classes[i])
                list = [j for j=start+1:finish]
            end
            push!(lists,list)
        end
        =#
    end
    permutations = get_all_permutations(lists)
    perm2 = []
    for permutation in permutations
        result = gwmt.permute_multiplication_table(permutation, multable)
        if all(result .== std_multable)
            perm2 = permutation
            break
        end
    end
    if perm2 == [] error("Error occurs at get_isomorphism_to_standard_point_group: can not find isomorphism between sorted group and standard point group") end

    final_perm = combine_permutations(perm1, perm2)
    multable = get_multiplication_table(group)
    result = gwmt.permute_multiplication_table(final_perm, multable)
    if !all(result .== std_multable) error("Error occurs at get_isomorphism_to_standard_point_group: find isorphism correspondence but wrong permutation") end

    return final_perm
end

"""
get_multiplication_table(pgnum::Int64) -> Matrix{Int64}
get_multiplication_table(encodeds::Vector{pg_encodetype}) -> Matrix{Int64}
get_multiplication_table(syms::Vector{pg_symtype}) -> Matrix{Int64}
"""
function get_multiplication_table(pgnum::Int64)::Matrix{Int64}
    syms = get_full_symmetry_operations(pgnum, encode=false)
    return get_multiplication_table(syms)
end
function get_multiplication_table(encodeds::Vector{pg_encodetype})::Matrix{Int64}
    syms = [pg_decode_symmetry(encoded) for encoded in encodeds]
    return get_multiplication_table(syms)
end
function get_multiplication_table(syms::Vector{pg_symtype})::Matrix{Int64}
    if !isa_group(syms) error("Error occurs at get_multiplication_table: symmetry operations don't form a group") end
    order = length(syms)
    mul_table = zeros(Int64, order, order)
    for i = 1:order
        for j = 1:order
            nsym = composite_symmetries(syms[i], syms[j], encode=false)
            idx = findall(sym->all(sym.==nsym), syms)[1]
            mul_table[i, j] = idx
        end
    end
    return mul_table
end

"""
conjugate_group(pgnum::Int64, conjg_sym_encoded::pg_encodetype; encode::Bool) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
conjugate_group(group::Vector{pg_encodetype}, conjg_sym_encoded::pg_encodetype; encode::Bool) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
conjugate_group(group::Vector{pg_symtype}, conjg_sym::pg_symtype; encode::Bool) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
"""
function conjugate_group(pgnum::Int64, conjg_sym_encoded::pg_encodetype; encode::Bool)::Union{Vector{pg_encodetype}, Vector{pg_symtype}}
    group = get_full_symmetry_operations(pgnum)
    conjg_sym = pg_decode_symmetry(conjg_sym_encoded)
    return conjugate_group(group, conjg_sym, encode=encode)
end

function conjugate_group(group::Vector{pg_encodetype}, conjg_sym_encoded::pg_encodetype; encode::Bool)::Union{Vector{pg_encodetype}, Vector{pg_symtype}}
    group = [pg_decode_symmetry(sym) for sym in group]
    conjg_sym = pg_decode_symmetry(conjg_sym_encoded)
    return conjugate_group(group, conjg_sym, encode=encode)
end

function conjugate_group(group::Vector{pg_symtype}, conjg_sym::pg_symtype; encode::Bool)::Union{Vector{pg_encodetype}, Vector{pg_symtype}}
    if encode
        conjg_group = Vector{pg_encodetype}()
    else
        conjg_group = Vector{pg_symtype}()
    end
    for ele in group
        push!(conjg_group, conjugate_symmetry(ele, conjg_sym, encode=encode))
    end
    return conjg_group
end

"""
get_class_of_point_group(syms::Union{Vector{pg_encodetype}, Vector{pg_symtype}}) -> Int64
get_class_of_point_group(syms::Union{pg_encodetype, pg_symtype}...) -> Int64
"""
function get_class_of_point_group(syms::Union{Vector{pg_encodetype}, Vector{pg_symtype}})::Int64
    return get_class_of_point_group(syms...)
end
function get_class_of_point_group(syms::Union{pg_encodetype, pg_symtype}...)::Int64
    iso_type_table = get_isometry_type_table(syms...)
    for pgnum in eachindex(pg_data)
        if all(iso_type_table .== pg_data[pgnum].iso_type_table) return pgnum end
    end
    error("Error occurs at get_class_of_point_group: can not find class of given point group")
end

"""
get_notation_of_point_group(syms::Union{Vector{pg_encodetype}, Vector{pg_symtype}})
get_notation_of_point_group(idx::Int64)
"""
function get_notation_of_point_group(syms::Union{Vector{pg_encodetype}, Vector{pg_symtype}})
    idx = get_class_of_point_group(syms)
    return get_notation_of_point_group(idx)
end
function get_notation_of_point_group(idx::Int64)
    symbol = pg_data[idx].HM_symbol
    return symbol
end


"""
get_isometry_type_table(syms::Union{Vector{pg_encodetype}, Vector{pg_symtype}}) -> Vector{Int64}
get_isometry_type_table(encodeds::pg_encodetype...) -> Vector{Int64}
get_isometry_type_table(syms::pg_symtype...) -> Vector{Int64}
isometry_type_table
index      1  2  3  4  5  6  7  8  9 10
operation -6 -4 -3 -2 -1  1  2  3  4  6
"""
function get_isometry_type_table(syms::Union{Vector{pg_encodetype}, Vector{pg_symtype}})::Vector{Int64}
    return get_isometry_type_table(syms...)
end
function get_isometry_type_table(encodeds::pg_encodetype...)::Vector{Int64}
    syms = [ pg_decode_symmetry(encoded) for encoded in encodeds]
    return get_isometry_type_table(syms...)
end
function get_isometry_type_table(syms::pg_symtype...)::Vector{Int64}
    error = zeros(Int64,10)
    table = zeros(Int64,10)
    for sym in syms
        iso_type = get_isometry_type(sym)
        if iso_type == 0
            return error
        else
            table[iso_type] += 1
        end
    end
    return table
end

"""
is_group_abelian(pgnum::Int64) -> Bool
check if a group is abelian
"""
function is_group_abelian(pgnum::Int64)::Bool
    encodeds = get_full_symmetry_operations(pgnum, encode=true)
    for i in eachindex(encodeds)
        if i == 1 continue end
        if !all([is_sym1_commuting_with_sym2(sym1=encodeds[i], sym2=encodeds[j]) for j = 1:i-1]) return false end
    end
    return true
end

"""
isa_set(list::Union{Vector{pg_encodetype}, Vector{pg_symtype}}) -> Bool
isa_set(encodeds::pg_encodetype...) -> Bool
isa_set(list::pg_symtype...) -> Bool
#check if a list forms a set,i.e.,each element only appears once
"""
function isa_set(list::Union{Vector{pg_encodetype}, Vector{pg_symtype}})::Bool
    return isa_set(list...)
end
function isa_set(encodeds::pg_encodetype...)::Bool
    list = [pg_decode_symmetry(encoded) for encoded in encodeds]
    return isa_set(list)
end
function isa_set(list::pg_symtype...)::Bool
    repeated = []
    for i in eachindex(list)
        for j in eachindex(list)
            if i == j continue end
            if list[i] != list[j]
                continue
            else
                if (i, j) in repeated || (j, i) in repeated continue end
                push!(repeated, (i, j))
            end
        end
    end
    if repeated == []
        return true
    else
        return false
    end
end

"""
isa_group(list::Union{Vector{pg_encodetype}, Vector{pg_symtype}}) -> Bool
isa_group(list::pg_encodetype...) -> Bool
isa_group(list::pg_symtype...) -> Bool
check if a set forms a pointgroup
"""
function isa_group(list::Union{Vector{pg_encodetype}, Vector{pg_symtype}})::Bool
    return isa_group(list...)
end
function isa_group(list::pg_encodetype...)::Bool
    list = [pg_decode_symmetry(ele) for ele in list]
    return isa_group(list...)
end
function isa_group(list::pg_symtype...)::Bool
    if !isa_set(list...)
        println("input has repeated matrix")
        return false
    end
    if !([1 0 0;0 1 0;0 0 1] in list) return false end
    for ele1 in list
        if !(Int64.(inv(ele1)) in list) return false end
        for ele2 in list
            if !(composite_symmetries(ele1, ele2, encode=false) in list) return false end
        end
    end
    return true
end

"""
get_full_symmetry_operations(pgnum::Int64; encode::Bool) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
#given a pointgroup number, get all its symmetry operations from PointGroupDataBase
"""
function get_full_symmetry_operations(pgnum::Int64; encode::Bool)::Union{Vector{pg_encodetype}, Vector{pg_symtype}}
    if pgnum <1 || pgnum >32 error("Error occurs at get_character_table: point group number should be in range of [1, 32]") end
    encodeds = [pg_symmetry_operations[i] for i = pg_index[pgnum]:pg_index[pgnum]+pg_data[pgnum].order-1]
    if encode
        return encodeds
    else
        syms = [pg_decode_symmetry(encoded) for encoded in encodeds]
        return syms
    end
end

"""
generate_full_symmetry_operations(syms::Union{Vector{pg_encodetype}, Vector{pg_symtype}}; encode::Bool, sort::Bool=false) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
generate_full_symmetry_operations(encodeds::pg_encodetype...; encode::Bool, sort::Bool=false) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
generate_full_symmetry_operations(syms::pg_symtype...; encode::Bool, sort::Bool=false) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
"""
function generate_full_symmetry_operations(syms::Union{Vector{pg_encodetype}, Vector{pg_symtype}}; encode::Bool, sort::Bool=false)::Union{Vector{pg_encodetype}, Vector{pg_symtype}}
    return generate_full_symmetry_operations(syms..., encode=encode, sort=sort)
end
function generate_full_symmetry_operations(encodeds::pg_encodetype...; encode::Bool, sort::Bool=false)::Union{Vector{pg_encodetype}, Vector{pg_symtype}}
    syms = [pg_decode_symmetry(encoded) for encoded in encodeds]
    return generate_full_symmetry_operations(syms..., encode=encode, sort=sort)
end
function generate_full_symmetry_operations(syms::pg_symtype...; encode::Bool, sort::Bool=false)::Union{Vector{pg_encodetype}, Vector{pg_symtype}}
    full_syms = push!(Vector{pg_symtype}(),[1 0 0;0 1 0;0 0 1])
    for sym in syms
        new_sym_list = [composite_symmetries(sym, order = i, encode=false) for i = 1:get_symmetry_order(sym)-1]
        while new_sym_list != []
            new_sym = new_sym_list[1]
            if !(new_sym in full_syms)
                push!(full_syms, new_sym)
                for fsym in full_syms
                    new_sym1 = composite_symmetries(fsym, new_sym, encode=false)
                    if !(new_sym1 ∈ full_syms) && !(new_sym1 ∈ new_sym_list) push!(new_sym_list, new_sym1) end
                    new_sym2 = composite_symmetries(new_sym, fsym, encode=false)
                    if !(new_sym2 ∈ full_syms) && !(new_sym2 ∈ new_sym_list) push!(new_sym_list, new_sym2) end
                end
            end
            new_sym_list = new_sym_list[2:end]
        end
    end

    #sort symmetry operations
    if sort
        conjg_classes = get_full_conjugacy_classes(full_syms, encode=encode)
        if encode
            full_syms = Vector{pg_encodetype}()
        else
            full_syms = Vector{pg_symtype}()
        end
        for class in conjg_classes
            append!(full_syms, class)
        end
        return full_syms
    elseif encode
        full_syms = [pg_encode_symmetry(sym) for sym in full_syms]
        return full_syms
    end

end

"""
get_full_subgroups(pgnum::Int64; encode::Bool) -> Union{Vector{Vector{pg_encodetype}}, Vector{Vector{pg_symtype}}}
#get all full subgroups of given group from PointGroupDataBase
"""
function get_full_subgroups(pgnum::Int64; encode::Bool)::Union{Vector{Vector{pg_encodetype}}, Vector{Vector{pg_symtype}}}
    subgroups = pg_full_subgroups[pgnum]
    if encode
        subgroups = [collect(subgroup) for subgroup in subgroups]
    else
        subgroups = [[pg_decode_symmetry(sym) for sym in subgroup] for subgroup in subgroups]
    end
    return subgroups
end

"""
get_certain_subgroup(pgnum::Int64, subindx::Int64; encode::Bool) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
#get a certain subgroup of given group from PointGroupDataBase
"""
function get_certain_subgroup(pgnum::Int64, subidx::Int64; encode::Bool)::Union{Vector{pg_encodetype}, Vector{pg_symtype}}
    if pgnum <1 || pgnum >32 error("Error occurs at get_character_table: point group number should be in range of [1, 32]") end
    subgroup = pg_full_subgroups[pgnum][subidx]
    if encode
        subgroup = collect(subgroup)
    else
        subgroup = [pg_decode_symmetry(sym) for sym in subgroup]
    end
    return subgroup
end

end#(module GroupInformation)

##################################################GroupRepresentation#################################################
######################################################################################################################
module GroupRepresentation
using LinearAlgebra
using ...MathTool
using ..PointGroupDataBase
using ..PointSymmetry
using ..GroupInformation

export get_irreps
export get_character_table
export get_regular_representation

#=
get_irreps(pgnum::Int64)::Matrix{Matrix{ComplexF64}}
get_irreps(group::Union{Vector{pg_encodetype}, Vector{pg_symtype}})::Matrix{Matrix{ComplexF64}}
get_character_table(pgnum::Int64) -> Matrix{Any}
get_character_table(group::Union{Vector{pg_encodetype}, Vector{pg_symtype}}) -> Matrix{Any}
get_regular_representation(pgnum::Int64) -> Vector{Matrix{Int64}}
get_regular_representation(syms::Vector{pg_encodetype}) -> Vector{Matrix{Int64}}
=#

"""
get_irreps(pgnum::Int64) -> Matrix{Matrix{ComplexF64}}
get_irreps(group::Union{Vector{pg_encodetype}, Vector{pg_symtype}}) -> Matrix{Matrix{ComplexF64}}
"""
function get_irreps(pgnum::Int64)::Matrix{Matrix{ComplexF64}}
    if pgnum < 1 || pgnum > 32 error("Error occurs at get_character_table: point group number should be in range of [1, 32]") end
    data = pg_irreps[pgnum]
    irrepnum = length(data)
    symnum = length(data[1])
    irreps = Matrix{Matrix{ComplexF64}}([Matrix{ComplexF64}([data[i][j];;]) for i in 1:irrepnum, j in 1:symnum])
    return irreps
end
function get_irreps(group::Union{Vector{pg_encodetype}, Vector{pg_symtype}})::Matrix{Matrix{ComplexF64}}
    if !isa_group(group) error("Error occurs at get_irreps: symmetry operations do not form a group") end
    std_pgnum = get_class_of_point_group(group...)
    std_group = get_full_symmetry_operations(std_pgnum, encode=true)
    if is_group_same(group, std_group)
        permutation = get_permutation_to_standard_point_group(group, std_group=std_group)
    else
        permutation = get_isomorphism_to_standard_point_group(group, std_group=std_group)
    end
    std_irreps = get_irreps(std_pgnum)
    irreps = reshape([std_irreps[j, idx] for idx in permutation for j = 1:size(std_irreps, 1)], size(std_irreps))
    return irreps
end

"""
get_character_table(pgnum::Int64) -> Matrix{Any}
get_character_table(group::Union{Vector{pg_encodetype}, Vector{pg_symtype}}) -> Matrix{Any}
"""
function get_character_table(pgnum::Int64)::Matrix{Any}
    if pgnum <1 || pgnum >32 error("Error occurs at get_character_table: point group number should be in range of [1, 32]") end
    table = pg_character_table[pgnum]
    irrepnum = length(table)
    symnum = length(table[1])
    character_table = reshape([table[j][i] for i = 1:symnum for j = 1:irrepnum], irrepnum, symnum)
    return character_table
end
function get_character_table(group::Union{Vector{pg_encodetype}, Vector{pg_symtype}})::Matrix{Any}
    if !isa_group(group) error("Error occurs at get_character_table: symmetry operations do not form a group") end

    std_pgnum = get_class_of_point_group(group...)
    std_group = get_full_symmetry_operations(std_pgnum, encode=true)
    if is_group_same(group, std_group)
        permutation = get_permutation_to_standard_point_group(group, std_group=std_group)
    else
        permutation = get_isomorphism_to_standard_point_group(group, std_group=std_group)
    end
    std_table = get_character_table(std_pgnum)
    character_table = reshape([std_table[j, idx] for idx in permutation for j = 1:size(std_table, 1)], size(std_table))
    return character_table
end

"""
get_regular_representation(pgnum::Int64) -> Vector{Matrix{Int64}}
get_regular_representation(syms::Vector{pg_encodetype}) -> Vector{Matrix{Int64}}
"""
function get_regular_representation(pgnum::Int64)::Vector{Matrix{Int64}}
    if pgnum <1 || pgnum >32 error("Error occurs at get_character_table: point group number should be in range of [1, 32]") end
    syms = get_full_symmetry_operations(pgnum, encode=false)
    return get_regular_representation(syms)
end
function get_regular_representation(syms::Vector{pg_encodetype})::Vector{Matrix{Int64}}
    order = length(syms)
    mul_table = zeros(Int64, order, order)
    syms_inv_idx = []
    for sym in syms
        sym_inv = get_inverse_of_matrix(pg_decode_symmetry(sym))
        if !(sym_inv in syms) error("Error:inverse of symmetry not in group") end
        idx = findall(x->x==sym_inv, syms)[1]
        push!(syms_inv_idx, idx)
    end
    for i = 1:order
        for j = 1:order
            nsym = composite_symmetries(syms[syms_inv_idx[i]], syms[j], encode=true)
            idx = findall(sym -> sym==nsym, syms)[1]
            mul_table[i, j] = idx
        end
    end
    reg_rep = Vector{Matrix{Int64}}()
    for idx = 1:order
        rep = reshape([i==idx ? 1 : 0 for i in mul_table], order, order)
        push!(reg_rep, rep)
    end
    return reg_rep
end

end#(module GroupRepresentation)

########################################################kSpace########################################################
######################################################################################################################
module kSpace
import ...kSpace as ksp
using LinearAlgebra
using ...MathTool
using ..PointGroupDataBase
using ..PointSymmetry
using ..GroupInformation

export get_high_symmetry_klist_of_group
export get_k_invariant_under_group
export get_k_invariant_under_symmetry_operation
export get_symmetry_in_kspace
export get_little_group_of_k
export is_k_invariant_under_group
export is_k_invariant_under_symmetry_operation
export is_k_equivalent_under_group

ktype = Matrix{Float64}
ksymtype = Matrix{Float64}

#=
get_symmetry_in_kspace(rencoded::pg_encodetype) -> ksymtype
get_symmetry_in_kspace(rsym::pg_symtype) -> ksymtype
get_high_symmetry_klist_of_group(idx::Int64) -> Dict{Int64, Vector{Any}}
get_k_invariant_under_group(idx::Int64) -> Vector{ktype}}
get_k_invariant_under_group(rencodeds::Vector{pg_encodetype}) -> Vector{ktype}
get_k_invariant_under_group(rsyms::Vector{pg_symtype}) -> Vector{ktype}
get_k_invariant_under_symmetry_operation(rencoded::pg_encodetype) -> Vector{ktype}
get_k_invariant_under_symmetry_operation(rsym::pg_symtype) -> Vector{ktype}
get_little_group_of_k(pg_num::Int64; k::ktype, encode::Bool) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
get_little_group_of_k(rencodeds::Vector{pg_encodetype}; k::ktype, encode::Bool) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
get_little_group_of_k(rsyms::Vector{pg_symtype}; k::ktype, encode::Bool) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
is_k_invariant_under_group(k::ktype, rencodeds::Vector{pg_encodetype}) -> Bool
is_k_invariant_under_group(k::ktype, rsyms::Vector{pg_symtype}) -> Bool
is_k_invariant_under_symmetry_operation(k::ktype, rencoded::pg_encodetype) -> Bool
is_k_invariant_under_symmetry_operation(k::ktype, rsym::pg_symtype) -> Bool
is_k_equivalent_under_group(k1::ktype, k2::ktype, rencodeds::Vector{pg_encodetype}) -> Bool
is_k_equivalent_under_group(k1::ktype, k2::ktype, rsyms::Vector{pg_symtype}) -> Bool
=#

"""
get_symmetry_in_kspace(rencoded::pg_encodetype) -> ksymtype
get_symmetry_in_kspace(rsym::pg_symtype) -> ksymtype
"""
function get_symmetry_in_kspace(rencoded::pg_encodetype)::ksymtype
    rsym = pg_decode_symmetry(rencoded)
    return get_symmetry_in_kspace(rsym)
end
function get_symmetry_in_kspace(rsym::pg_symtype)::ksymtype
    ksym = transpose(inv(rsym))
    return ksym
end

"""
get_high_symmetry_klist_of_group(idx::Int64) -> Dict{Int64, Vector{Any}}
"""
function get_high_symmetry_klist_of_group(idx::Int64)::Dict{Int64, Vector{Any}}
    if !(idx in 1:32) error("Error occurs at get_high_symmetry_klist_of_group: point group index out of range [1, 32]") end
    group = get_full_symmetry_operations(idx, encode=true)
    subgroups = get_full_subgroups(idx, encode=true)
    klist = []
    for i in eachindex(subgroups)
        rencodeds = subgroups[i]
        for k in get_k_invariant_under_group(rencodeds)
            if !ksp.is_k_in_klist(k, klist) push!(klist, k) end
        end
    end
    highsymk = Dict(0=>[], 1=>[], 2=>[], 3=>[])
    for k in klist
        n = ksp.get_k_geometric_type(k)
        for i in eachindex(highsymk[n])
            ele = highsymk[n][i]
            if is_k_equivalent_under_group(k, ele[1], group)
                push!(highsymk[n][i], k)
                @goto next_k
            end
        end
        push!(highsymk[n],[k])
        @label next_k
    end

    return highsymk
end

"""
get_k_invariant_under_group(idx::Int64) -> Vector{ktype}}
get_k_invariant_under_group(rencodeds::Vector{pg_encodetype}) -> Vector{ktype}
get_k_invariant_under_group(rsyms::Vector{pg_symtype}) -> Vector{ktype}
"""
function get_k_invariant_under_group(idx::Int64)::Vector{ktype}
    rsyms = get_full_symmetry_operations(idx, encode=false)
    return get_k_invariant_under_group(rsyms)
end
function get_k_invariant_under_group(rencodeds::Vector{pg_encodetype})::Vector{ktype}
    rsyms = [pg_decode_symmetry(rencoded) for rencoded in rencodeds]
    return get_k_invariant_under_group(rsyms)
end
function get_k_invariant_under_group(rsyms::Vector{pg_symtype})::Vector{ktype}
    ksyms = [get_symmetry_in_kspace(rsym) for rsym in rsyms]
    k_list = ksp.get_k_invariant_under_group(ksyms)
    for k in k_list
        encodeds = get_little_group_of_k(rsyms, k=k, encode=true)
        if length(encodeds) != length(rsyms) || !is_k_invariant_under_group(k, rsyms)
            error("Error occurs at get_k_invariant_under_group: get wrong k_list")
        end
    end
    return k_list
end

"""
get_k_invariant_under_symmetry_operation(rencoded::pg_encodetype) -> Vector{ktype}
get_k_invariant_under_symmetry_operation(rsym::pg_symtype) -> Vector{ktype}
"""
function get_k_invariant_under_symmetry_operation(rencoded::pg_encodetype)::Vector{ktype}
    rsym = pg_decode_symmetry(rencoded)
    return get_k_invariant_under_symmetry_operation(rsym)
end
function get_k_invariant_under_symmetry_operation(rsym::pg_symtype)::Vector{ktype}
    ksym = get_symmetry_in_kspace(rsym)
    return ksp.get_k_invariant_under_symmetry_operation(ksym)
end

"""
get_little_group_of_k(pg_num::Int64; k::ktype, encode::Bool) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
get_little_group_of_k(rencodeds::Vector{pg_encodetype}; k::ktype, encode::Bool) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
get_little_group_of_k(rsyms::Vector{pg_symtype}; k::ktype, encode::Bool) -> Union{Vector{pg_encodetype}, Vector{pg_symtype}}
#given point group and k, find little group at k
"""
function get_little_group_of_k(pg_num::Int64; k::ktype, encode::Bool)::Union{Vector{pg_encodetype}, Vector{pg_symtype}}
    rencodeds = get_full_symmetry_operations(pg_num, encode=true)
    return get_little_group_of_k(rencodeds, k=k, encode=encode)
end
function get_little_group_of_k(rencodeds::Vector{pg_encodetype}; k::ktype, encode::Bool)::Union{Vector{pg_encodetype}, Vector{pg_symtype}}
    rsyms = [pg_decode_symmetry(rencoded) for rencoded in rencodeds]
    return get_little_group_of_k(rsyms, k=k, encode=encode)
end
function get_little_group_of_k(rsyms::Vector{pg_symtype}; k::ktype, encode::Bool)::Union{Vector{pg_encodetype}, Vector{pg_symtype}}
    little_group = Vector{pg_symtype}()
    for idx in eachindex(rsyms)
        ksym = get_symmetry_in_kspace(rsyms[idx])
        #println(ksym)
        if ksp.is_k_invariant_under_symmetry_operation(k, ksym) push!(little_group, rsyms[idx]) end
    end
    if encode little_group = [pg_encode_symmetry(rsym) for rsym in little_group] end
    return little_group
end

"""
is_k_invariant_under_group(k::ktype, rencodeds::Vector{pg_encodetype}) -> Bool
is_k_invariant_under_group(k::ktype, rsyms::Vector{pg_symtype}) -> Bool
"""
function is_k_invariant_under_group(k::ktype, rencodeds::Vector{pg_encodetype})::Bool
    rsyms = [pg_decode_symmetry(rencoded) for rencoded in rencodeds]
    return is_k_invariant_under_group(k, rsyms)
end
function is_k_invariant_under_group(k::ktype, rsyms::Vector{pg_symtype})::Bool
    ksyms = [get_symmetry_in_kspace(rsym) for rsym in rsyms]
    return ksp.is_k_invariant_under_group(k, ksyms)
end

"""
is_k_invariant_under_symmetry_operation(k::ktype, rencoded::pg_encodetype) -> Bool
is_k_invariant_under_symmetry_operation(k::ktype, rsym::pg_symtype) -> Bool
"""
function is_k_invariant_under_symmetry_operation(k::ktype, rencoded::pg_encodetype)::Bool
    rsym = pg_decode_symmetry(rencoded)
    return is_k_invariant_under_symmetry_operation(k, rsym)
end
function is_k_invariant_under_symmetry_operation(k::ktype, rsym::pg_symtype)::Bool
    ksym = get_symmetry_in_kspace(rsym)
    return ksp.is_k_invariant_under_symmetry_operation(k, ksym)
end

"""
is_k_equivalent_under_group(k1::ktype, k2::ktype, rencodeds::Vector{pg_encodetype}) -> Bool
is_k_equivalent_under_group(k1::ktype, k2::ktype, rsyms::Vector{pg_symtype}) -> Bool
"""
function is_k_equivalent_under_group(k1::ktype, k2::ktype, rencodeds::Vector{pg_encodetype})::Bool
    rsyms = [pg_decode_symmetry(rencoded) for rencoded in rencodeds]
    return is_k_equivalent_under_group(k1, k2, rsyms)
end
function is_k_equivalent_under_group(k1::ktype, k2::ktype, rsyms::Vector{pg_symtype})::Bool
    ksyms = [get_symmetry_in_kspace(rsym) for rsym in rsyms]
    return ksp.is_k_equivalent_under_group(k1, k2, ksyms)
end

end#(module kSpace)

using .PointGroupDataBase
using .PointSymmetry
using .GroupInformation
using .GroupRepresentation
using .kSpace
end#(module PointGroup)
