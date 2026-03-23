try
    typeof(GroupWithMultiplicationTable)
catch UndefVarError
    include("../representation-tool/groupwithmultiplicationtable.jl")
end
try
    typeof(Tool)
catch UndefVarError
    include("../math-tool/tool.jl")
end
######################################################################################################################
####################################################AbstractGroup#####################################################
######################################################################################################################
module AbstractGroup
include("../database/abg_database.jl")

####################################################SpaceSymmetry#####################################################
######################################################################################################################
module AbstractSymmetry
using ..AbstractGroupDataBase
import ...Tool

export is_relation_commutative
export is_symmetry_contained
export is_symmetry_equivalent
export get_order_of_generator
export get_order_of_symmetry
export composite_symmetries
export sort_symmetry_in_order_of_generators
export exchange_first_noncommutative_neighbor_generators
export exchange_commutative_neighbor_generators
export remove_redundant_generator
export is_symmetry_in_order_generators
export conjugate_symmetry
export get_inverse_of_symmetry

"""
is_relation_commutative(sym1::Vector{Int64}, sym2::Vector{Int64}) -> Bool
"""
function is_relation_commutative(sym1::Vector{Int64}, sym2::Vector{Int64})::Bool
    if length(sym1) != length(sym2) return false end
    len = length(sym1)
    for idx in 1:len-1
        len1 = idx
        len2 = len-len1
        if sym1[1:len1] == sym2[len-len1+1:len] && sym1[len-len2+1:len] == sym2[1:len2] return true end
    end
    return false
end

"""
is_symmetry_contained(symmetry::Vector{Int64}, symmetries::Vector{Vector{Int64}}; relns::abg_relns) -> Bool
"""
function is_symmetry_contained(symmetry::Vector{Int64}, symmetries::Vector{Vector{Int64}}; relns::abg_relns)::Bool
    for sym in symmetries
        if is_symmetry_equivalent(sym, symmetry, relns=relns) return true end
    end
    return false
end

"""
is_symmetry_equivalent(symmetry1::Vector{Int64}, symmetry2::Vector{Int64}; relns::abg_relns) -> Bool
"""
function is_symmetry_equivalent(symmetry1::Vector{Int64}, symmetry2::Vector{Int64}; relns::abg_relns)::Bool
    sym1 = deepcopy(symmetry1)
    sym2 = deepcopy(symmetry2)
    #sym1 = remove_redundant_generator(sym1, relns=relns)
    #sym2 = remove_redundant_generator(sym2, relns=relns)
    if sym1 == sym2 return true end
    return false
    order1 = get_order_of_symmetry(sym1, relns=relns)
    order2 = get_order_of_symmetry(sym2, relns=relns)
    if order1 != order2 return false end
    for i in 1:order1-1
        j = order1-i
        isym1 = composite_symmetries(sym1, order=i, relns=relns)
        jsym2 = composite_symmetries(sym2, order=j, relns=relns)
        if composite_symmetries(isym1, jsym2, relns=relns) != [1] return false end
    end
    return true
end

"""
get_order_of_generator(generator::Int64; relns::abg_relns) -> Int64
"""
function get_order_of_generator(generator::Int64; relns::abg_relns)::Int64
    gen = deepcopy(generator)
    for pair in relns[1]
        if pair[2] == [] && gen in pair[1] return length(pair[1]) end
    end
    error("can not get order of generator")
end

"""
get_order_of_symmetry(symmetry::Vector{Int64}; relns::abg_relns) -> Int64
"""
function get_order_of_symmetry(symmetry::Vector{Int64}; relns::abg_relns)::Int64
    sym = deepcopy(symmetry)
    for i in 1:20
        if composite_symmetries(sym,order=i,relns=relns) == [1] return i end
    end
    error("can not get order of symmetry")
end

"""
composite_symmetries(sym::Vector{Int64}; order::Int64=1, relns::abg_relns) -> Vector{Int64}
composite_symmetries(syms::Vector{Int64}...; relns::abg_relns) -> Vector{Int64}
"""
function composite_symmetries(sym::Vector{Int64}; order::Int64=1, relns::abg_relns)::Vector{Int64}
    if order == 0
        return [1]
    elseif order == 1
        return sort_symmetry_in_order_of_generators(sym, relns=relns)
    else
        return composite_symmetries([sym for i = 1:order]..., relns=relns)
    end
end
function composite_symmetries(syms::Vector{Int64}...; relns::abg_relns)::Vector{Int64}
    if length(syms) == 1
        return sort_symmetry_in_order_of_generators(syms[1], relns=relns)
    elseif length(syms) == 2
        sym = [syms[1]; syms[2]]
        return sort_symmetry_in_order_of_generators(sym, relns=relns)
    else
        return composite_symmetries(composite_symmetries(syms[1], syms[2], relns=relns), syms[3:end]..., relns=relns)
    end
end

"""
sort_symmetry_in_order_of_generators(symmetry::Vector{Int64}; relns::abg_relns) -> Vector{Int64}
"""
function sort_symmetry_in_order_of_generators(symmetry::Vector{Int64}; relns::abg_relns)::Vector{Int64}
    sym = deepcopy(symmetry)
    sym = remove_redundant_generator(sym, relns=relns)
    step = 1
    while !is_symmetry_in_order_generators(sym, relns=relns)
        if step > 1e5 error("step larger that 1e5 in sort_symmetry_in_order_of_generators") end
        sym = exchange_commutative_neighbor_generators(sym, relns=relns)
        sym = exchange_first_noncommutative_neighbor_generators(sym, relns=relns)
        step+=1
    end
    if length(sym) == 0 sym = [1] end
    return sym
end

"""
exchange_first_noncommutative_neighbor_generators(symmetry::Vector{Int64}; relns::abg_relns) -> Vector{Int64}
"""
function exchange_first_noncommutative_neighbor_generators(symmetry::Vector{Int64}; relns::abg_relns)::Vector{Int64}
    sym = deepcopy(symmetry)
    firstidxs = nothing
    reln = nothing
    for pair in relns[4]
        idxs = Tool.find_range_of_subarray(pair[1], sym)
        if length(idxs) == 0 continue end
        if firstidxs == nothing
            firstidxs = idxs
            reln = pair[2]
        elseif idxs[1] < firstidxs[1]
            firstidxs = idxs
            reln = pair[2]
        end
    end
    if firstidxs == nothing return sym end
    sym = [sym[1:firstidxs[end]]; reln; sym[firstidxs[end]+1:end]]
    deleteat!(sym, firstidxs)
    return remove_redundant_generator(sym, relns=relns)
end

"""
exchange_commutative_neighbor_generators(symmetry::Vector{Int64}; relns::abg_relns) -> Vector{Int64}
"""
function exchange_commutative_neighbor_generators(symmetry::Vector{Int64}; relns::abg_relns)::Vector{Int64}
    sym = deepcopy(symmetry)
    for pair in relns[3]
        idxs = Tool.find_range_of_subarray(pair[1], sym)
        step = 1
        while length(idxs) != 0
            if step > 1e5 error("step larger that 1e5 in exchange_commutative_neighbor_generators") end
            sym = [sym[1:idxs[end]]; pair[2]; sym[idxs[end]+1:end]]
            deleteat!(sym, idxs)
            idxs = Tool.find_range_of_subarray(pair[1], sym)
            step+=1
        end
    end
    return remove_redundant_generator(sym, relns=relns)
end

"""
remove_redundant_generator(symmetry::Vector{Int64}; relns::abg_relns) -> Vector{Int64}
"""
function remove_redundant_generator(symmetry::Vector{Int64}; relns::abg_relns)::Vector{Int64}
    sym = deepcopy(symmetry)
    step = 1
    @label remove
    if step > 1e5 error("step larger that 1e5 in remove_redundant_generator") end
    for pair in relns[1]
        idxs = Tool.find_range_of_subarray(pair[1], sym)
        step1 = 1
        while length(idxs) != 0
            if step1 > 1e5 error("step1 larger that 1e5 in remove_redundant_generatorrs") end
            sym = [sym[1:idxs[end]]; pair[2]; sym[idxs[end]+1:end]]
            deleteat!(sym, idxs)
            idxs = Tool.find_range_of_subarray(pair[1], sym)
            step1+=1
        end
    end
    for pair in relns[2]
        idxs = Tool.find_range_of_subarray(pair[1], sym)
        step2=1
        while length(idxs) != 0
            if step2 > 1e5 error("step2 larger that 1e5 in remove_redundant_generator") end
            sym = [sym[1:idxs[end]]; pair[2]; sym[idxs[end]+1:end]]
            deleteat!(sym, idxs)
            idxs = Tool.find_range_of_subarray(pair[1], sym)
            step2+=1
        end
    end
    #=check=#
    for pair in relns[1]
        idxs = Tool.find_range_of_subarray(pair[1], sym)
        if length(idxs) == 0  continue end
        step+=1
        @goto remove
    end
    for pair in relns[2]
        idxs = Tool.find_range_of_subarray(pair[1], sym)
        if length(idxs) == 0 continue end
        step+=1
        @goto remove
    end
    if length(sym) == 0 sym = [1] end
    return sym
end

"""
is_symmetry_in_order_generators(symmetry::Vector{Int64}; relns::abg_relns) -> Bool
"""
function is_symmetry_in_order_generators(symmetry::Vector{Int64}; relns::abg_relns)::Bool
    sym = deepcopy(symmetry)
    for idx in eachindex(sym)
        for i in 1:idx-1
            if sym[i] > sym[idx] return false end
        end
    end
    return true
end

"""
conjugate_symmetry(sym::Vector{Int64}, conjg_sym::Vector{Int64}; relns::abg_relns) -> Vector{Int64}
"""
function conjugate_symmetry(sym::Vector{Int64}, conjg_sym::Vector{Int64}; relns::abg_relns)::Vector{Int64}
    return composite_symmetries(get_inverse_of_symmetry(conjg_sym, relns=relns), sym, conjg_sym, relns=relns)
end

"""
get_inverse_of_symmetry(symmetry::Vector{Int64}; relns::abg_relns) -> Vector{Int64}
"""
function get_inverse_of_symmetry(symmetry::Vector{Int64}; relns::abg_relns)::Vector{Int64}
    sym = deepcopy(symmetry)
    order = get_order_of_symmetry(sym, relns=relns)
    return composite_symmetries(sym,order=order-1,relns=relns)
end

end #module AbstractSymmetry
####################################################SpaceSymmetry#####################################################
######################################################################################################################
module GroupInformation
using ..AbstractGroupDataBase
using ..AbstractSymmetry
import ...GroupWithMultiplicationTable as gwmt

export get_full_conjugacy_classes
export get_multiplication_table
export isa_group

"""
get_multiplication_table(syms::Vector{Vector{Int64}}; relns::abg_relns) -> Matrix{Int64}
"""
function get_multiplication_table(syms::Vector{Vector{Int64}}; relns::abg_relns)::Matrix{Int64}
    #if !isa_group(syms, relns=relns, gens=gens) error("Error occurs at get_multiplication_table: symmetry operations don't form a group") end
    order = length(syms)
    mul_table = zeros(Int64, order, order)
    for i = 1:order
        for j = 1:order
            if false
                println(i, " ",j)
                println("symi=$(syms[i]),symj=$(syms[j])")
                nsym = composite_symmetries(syms[i], syms[j], relns=relns)
                println("nsym=$nsym")
                println(findall(sym->is_symmetry_equivalent(sym, nsym, relns=relns), syms))
                if length(findall(sym->is_symmetry_equivalent(sym, nsym, relns=relns), syms)) != 1
                    error("can not build multiplication table")
                end
                idx = findall(sym->is_symmetry_equivalent(sym, nsym, relns=relns), syms)[1]
                println("idx=$idx")
                mul_table[i, j] = idx
            else
                nsym = composite_symmetries(syms[i], syms[j], relns=relns)
                if length(findall(sym->is_symmetry_equivalent(sym, nsym, relns=relns), syms)) != 1
                    error("can not build multiplication table")
                end
                idx = findall(sym->is_symmetry_equivalent(sym, nsym, relns=relns), syms)[1]
                mul_table[i, j] = idx
            end
        end
    end
    return mul_table
end

#=
"""
get_full_conjugacy_classes(syms::Vector{abg_symtype}; relns, gens::Vector{abg_gentype}) -> Vector{Vector{abg_symtype}}
get_full_conjugacy_classes(syms::Vector{abg_symtype}; multable::Matrix{Int64}) -> Vector{Vector{abg_symtype}}
"""
function get_full_conjugacy_classes(syms::Vector{abg_symtype}; relns, gens::Vector{abg_gentype})::Vector{Vector{abg_symtype}}
    syms_remain = deepcopy(syms)
    conjg_classes = Vector{abg_symtype}[]
    while syms_remain != []
        class = abg_symtype[]
        push!(class, syms_remain[1])
        for sym in syms
            nsym = conjugate_symmetry(syms_remain[1], sym, relns=relns, gens=gens)
            if !(nsym in syms) error("Error occurs at get_conjugacy_classes_of_point_groups: group not closed") end
            if !(nsym in class) push!(class, nsym) end
        end
        push!(conjg_classes, class)
        syms_remain = [sym for sym in syms_remain if !(sym in class)]
    end
    return conjg_classes
end
function get_full_conjugacy_classes(syms::Vector{abg_symtype}; multable::Matrix{Int64})::Vector{Vector{abg_symtype}}
    classes_idx = gwmt.get_full_conjugacy_classes(multable=multable)
    classes = []
    for class_idx in classes_idx
        class = [syms[idx] for idx in class_idx]
        push!(classes, class)
    end
    return classes
end

"""
isa_group(full_syms::Vector{abg_symtype}; relns, gens::Vector{abg_symtype}) -> Bool
"""
function isa_group(full_syms::Vector{abg_symtype}; relns, gens::Vector{abg_gentype})::Bool
    #=check composition=#
    for sym1 in full_syms
        for sym2 in full_syms
            nsym = composite_symmetries(sym1, sym2, relns=relns, gens=gens)
            if !(nsym ∈ full_syms)
                println(sym1, sym2)
                println(nsym)
                return false
            end
        end
    end
    #=check inversion=#
    for sym1 in full_syms
        for sym2 in full_syms
            nsym = composite_symmetries(sym1, sym2, relns=relns, gens=gens)
            if nsym == full_syms[1] @goto next_sym1 end
        end
        return false
        @label next_sym1
    end
    return true
end
=#

end # module GroupInformation

using .AbstractGroupDataBase
using .AbstractSymmetry
using .GroupInformation
end #module AbstractGgroup
