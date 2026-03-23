try
    typeof(MathTool)
catch UndefVarError
    include("../math-tool/mathtool.jl")
end

module GroupWithMultiplicationTable
using LinearAlgebra
using Combinatorics
import ..MathTool

export get_isomorphism
export get_full_conjugacy_classes
export conjugate_symmetry
export composite_symmetries
export get_inverse_of_symmetry
export get_index_of_identity_element
export get_order_of_symmetry
export get_a_minimal_generating_set_of_group
export isa_normal_subgroup
export isa_subgroup
export generate_full_elements_in_order_of_generators
export generate_full_elements
export isa_group
export sort_symmetry_in_order_of_generators
export is_group_same
export get_generating_relations
export is_generating_relation_same
export permute_multiplication_table


"""
    permutation:P[i]
    syms1 = [1,2,3,...,i]
    syms2 = [,....,P[i],]
"""
function get_isomorphism(multable1::Matrix{Int64}, multable2::Matrix{Int64})
    genset1 = get_a_minimal_generating_set_of_group(multable=multable1)
    orderlist1 = [get_order_of_symmetry(gen, multable=multable1) for gen in genset1]
    syms1 = Vector{Int64}[]
    for idx in axes(multable1, 1)
        sym = sort_symmetry_in_order_of_generators(idx, gens=genset1, orderlist=orderlist1, multable=multable1)
        if sym == [] error(0) end
        push!(syms1, sym)
    end

    len = length(genset1)
    Eidx = get_index_of_identity_element(multable=multable2)
    idxs = [i for i in axes(multable2, 1) if i != Eidx]
    for perm in permutations(idxs, len)
        if length(Set(generate_full_elements_in_order_of_generators(perm, multable=multable2))) != size(multable2, 1)
            continue
        end
        orderlist2 = [get_order_of_symmetry(gen, multable=multable2) for gen in perm]
        if orderlist2 != orderlist1 continue end
        syms2 = Vector{Int64}[]
        for idx in axes(multable2, 1)
            sym2 = sort_symmetry_in_order_of_generators(idx, gens=perm, orderlist=orderlist2, multable=multable2)
            if sym2 == [] error(2) end
            push!(syms2, sym2)
        end
        if !is_group_same(syms1, syms2) continue end
        permutation = Vector{Int64}([findall(x->x==sym, syms2)[1] for sym in syms1])
        if all(multable2-permute_multiplication_table(permutation, multable1).==0)
            return permutation
        end
    end
    return Int64[]
end

"""
get_full_conjugacy_classes(multiplication_table::Matrix{Int64}) -> Vector{Vector{Int64}}
"""
function get_full_conjugacy_classes(;multable::Matrix{Int64})::Vector{Vector{Int64}}
    table = deepcopy(multable)
    idxs = [i for i in axes(table, 1)]
    idxs_remain = [i for i in axes(table, 1)]
    conjg_classes = Vector{Int64}[]
    while idxs_remain != []
        class = Int64[]
        push!(class, idxs_remain[1])
        for idx in idxs
            nidx = conjugate_symmetry(idxs_remain[1], idx, multable=multable)
            if !(nidx in idxs) error("Error occurs at get_conjugacy_classes: group not closed") end
            if !(nidx in class) push!(class, nidx) end
        end
        push!(conjg_classes, class)
        idxs_remain = [idx for idx in idxs_remain if !(idx in class)]
    end
    return conjg_classes
end

"""
conjugate_symmetry(idx::Int64, conjg_idx::Int64; multable::Matrix{Int64}) -> Int64
"""
function conjugate_symmetry(idx::Int64, conjg_idx::Int64; multable::Matrix{Int64})::Int64
    return composite_symmetries(get_inverse_of_symmetry(conjg_idx, multable=multable), idx, conjg_idx, multable=multable)
end

"""
composite_symmetries(idx::Int64; order::Int64, multable::Matrix{Int64}) -> Int64
composite_symmetries(idxs::Int64...; multable::Matrix{Int64}) -> Int64
"""
function composite_symmetries(idx::Int64; order::Int64=1, multable::Matrix{Int64})::Int64
    if order == 0
        return get_index_of_identity_element(multable=multable)
    elseif order == 1
        return idx
    else
        return composite_symmetries([idx for i in 1:order]..., multable=multable)
    end
end
function composite_symmetries(idxs::Int64...; multable::Matrix{Int64})::Int64
    if length(idxs) == 1
        return idxs[1]
    elseif length(idxs) == 2
        return multable[idxs[1], idxs[2]]
    else
        return composite_symmetries(composite_symmetries(idxs[1], idxs[2], multable=multable), idxs[3:end]..., multable=multable)
    end
end
function composite_symmetries(idxs::Vector{Int64}, list::Vector{Int64}; multable::Matrix{Int64})
    return composite_symmetries([composite_symmetries(idxs[i], order=list[i], multable=multable) for i in eachindex(list)]..., multable=multable)
end


"""
get_inverse_of_symmetry(idx::Int64; multable::Matrix{Int64}) -> Int64
"""
function get_inverse_of_symmetry(idx::Int64; multable::Matrix{Int64})::Int64
    table = deepcopy(multable)
    Eidx = get_index_of_identity_element(multable=table)
    inv_idx = 0
    for i in axes(table, 1)
        if table[i, idx] == Eidx
            inv_idx = i
            break
        end
    end
    if inv_idx == 0 error("can not find inverse of symmetry") end
    return inv_idx

end

"""
get_index_of_identity_element(;multable::Matrix{Int64}) -> Int64
"""
function get_index_of_identity_element(;multable::Matrix{Int64})::Int64
    Eidx = 0
    for i in axes(multable, 1)
       if multable[i, :] == [j for j in axes(multable, 2)]
            Eidx = i
            break
        end
    end
    if Eidx == 0 error("can not find identity symmetry") end
    return Eidx
end

"""
get_order_of_symmetry(idx::Int64; multable::Matrix{Int64}) -> Int64
"""
function get_order_of_symmetry(idx::Int64; multable::Matrix{Int64})::Int64
    Eidx = get_index_of_identity_element(multable=multable)
    order = 1
    while composite_symmetries(idx, order=order, multable=multable) != Eidx
        order+=1
        if order > 20 error("can not find order of symmetry") end
    end
    return order
end

"""
get_a_minimal_generating_set_of_group(;multable::Matrix{Int64}) -> Vector{Int64}
"""
function get_a_minimal_generating_set_of_group(;multable::Matrix{Int64})::Vector{Int64}
    len = floor(Int64, log(2, size(multable, 1))+1)
    Eidx = get_index_of_identity_element(multable=multable)
    idxs = [i for i in axes(multable, 1) if i != Eidx]
    if len == 1 return [1] end
    for i in 1:len
        for perm in permutations(idxs, i)
            if length(Set(generate_full_elements_in_order_of_generators(perm, multable=multable))) == size(multable, 1)
                return perm
            end
        end
    end
end

"""
isa_normal_subgroup(subgroup::Vector{Int64}, group::Vector{Int64})::Bool
"""
function isa_normal_subgroup(subgroup::Vector{Int64}, group::Vector{Int64}; multable::Matrix{Int64})::Bool
    for idx1 in subgroup
        for idx2 in group
            if !(conjugate_symmetry(idx1, idx2, multable=multable) in subgroup) return false end
        end
    end
    return true
end

"""
isa_subgroup(subgroup::Vector{Int64}, group::Vector{Int64}) -> Bool
"""
function isa_subgroup(subgroup::Vector{Int64}, group::Vector{Int64})::Bool
    for idx in subgroup
        if !(idx in group) return false end
    end
    return true
end

"""
generate_full_elements_in_order_of_generators(gens::Vector{Int64}; multable::Matrix{Int64}) -> Vector{Int64}
"""
function generate_full_elements_in_order_of_generators(gens::Vector{Int64}; multable::Matrix{Int64})::Vector{Int64}
    orderlist = [get_order_of_symmetry(gen, multable=multable) for gen in gens]
    syms=Int64[]
    for list in MathTool.generate_list(orderlist)
        nsym = composite_symmetries(gens, list, multable=multable)
        if !(nsym ∈ syms) push!(syms, nsym) end
    end
    return syms
end


"""
generate_full_elements(gens::gens::Vector{Int64}; multable::Matrix{Int64}) -> Vector{Int64}
"""
function generate_full_elements(gens::Vector{Int64}; multable::Matrix{Int64})::Vector{Int64}
    fullidxs = Int64[]
    nfullidxs = Int64[[1]; gens]
    while length(nfullidxs) != length(fullidxs)
        fullidxs = deepcopy(nfullidxs)
        nfullidxs = Int64[]
        for idxi in fullidxs, idxj in fullidxs
            nidx = composite_symmetries(idxi, idxj, multable=multable)
            if nidx in nfullidxs continue end
            push!(nfullidxs, nidx)
        end
    end
    if !isa_group(nfullidxs, multable=multable)
        println(nfullidxs)
        error("not a group")
    end
    return nfullidxs
    #=
    for gen in gens
        order = get_order_of_symmetry(gen, multable=multable)
        for i = 1:order-1
            sym_i = composite_symmetries(gen, order=i, multable=multable)
            if !(sym_i in full_elements) push!(full_elements, sym_i) end
            k = 1
            while k <= length(full_elements)
                nsym1 = composite_symmetries(full_elements[k], sym_i, multable=multable)
                if !(nsym1 in full_elements)
                    order1 = get_order_of_symmetry(nsym1, multable=multable)
                    for j = 1:order1-1
                        nsym1_j = composite_symmetries(nsym1, order=j, multable=multable)
                        if !(nsym1_j in full_elements) push!(full_elements, nsym1_j) end
                    end
                end
                nsym2 = composite_symmetries(sym_i, full_elements[k], multable=multable)
                if !(nsym2 in full_elements)
                    order2 = get_order_of_symmetry(nsym2, multable=multable)
                    for j = 1:order2-1
                        nsym2_j = composite_symmetries(nsym2, order=j, multable=multable)
                        if !(nsym2_j in full_elements) push!(full_elements, nsym2_j) end
                    end
                end
                k+=1
            end
        end
    end

    if !isa_group(full_elements, multable=multable)
        println(full_elements)
        error("not a group")
    end
    return full_elements
    =#
end


"""
isa_group(idxs::Vector{Int64}; multable::Matrix{Int64}) -> Bool
"""
function isa_group(idxs::Vector{Int64}; multable::Matrix{Int64})::Bool
    #=check composition=#
    for idx1 in idxs
        for idx2 in idxs
            nidx= composite_symmetries(idx1, idx2, multable=multable)
            if !(nidx ∈ idxs)
                return false
            end
        end
    end
    #=check inversion=#
    Eidx = get_index_of_identity_element(multable=multable)
    for idx1 in idxs
        for idx2 in idxs
            nidx = composite_symmetries(idx1, idx2, multable=multable)
            if nidx == Eidx @goto next_idx1 end
        end
        return false
        @label next_idx1
    end
    return true
end

"""
sort_symmetry_in_order_of_generators(idx::Int64; gens::Vector{Int64}, orderlist::Vector{Int64}, multable::Matrix{Int64}) -> Vector{Int64}
"""
function sort_symmetry_in_order_of_generators(idx::Int64; gens::Vector{Int64}, orderlist::Vector{Int64}, multable::Matrix{Int64})::Vector{Int64}
    #println("idx=",idx)
    for list in MathTool.generate_list(orderlist)
        #println(composite_symmetries(gens, list, multable=multable))
        if idx == composite_symmetries(gens, list, multable=multable)
            return list
        end
    end
    return []
    error("can not sort_symmetry_in_order_of_generators")
end

"""
is_group_same(idxs1::Vector{Vector{Int64}}, idxs2::Vector{Vector{Int64}}) -> Bool
"""
function is_group_same(idxs1::Vector{Vector{Int64}}, idxs2::Vector{Vector{Int64}})::Bool
    if length(idxs1) != length(idxs2) return false end
    for idx1 in idxs1
        if !(idx1 in idxs2) return false end
    end
    for idx2 in idxs2
        if !(idx2 in idxs1) return false end
    end
    return true
end

"""
get_generating_relations(generators::Vector{Int64}; multable::Matrix{Int64}) -> Dict{Vector{Int64}, Vector{Int64}}
"""
function get_generating_relations(generators::Vector{Int64}; multable::Matrix{Int64})::Dict{Vector{Int64}, Vector{Int64}}
    gens = deepcopy(generators)
    orderlist = [get_order_of_symmetry(gen, multable=multable) for gen in gens]
    revgens = reverse(gens)
    revorderlist = [get_order_of_symmetry(gen, multable=multable) for gen in revgens]
    relns = Dict{Vector{Int64}, Vector{Int64}}()
    for revlist in MathTool.generate_list(revorderlist)
        if sum(revlist) <= 1 continue end
        idx = composite_symmetries(revgens, revlist, multable=multable)
        relns[revlist] = sort_symmetry_in_order_of_generators(idx, gens=gens, orderlist=orderlist, multable=multable)
    end
    return relns
end

"""
is_generating_relation_same(relns1, relns2) -> Bool
"""
function is_generating_relation_same(relns1, relns2)::Bool
    for pair1 in relns1
        for pair2 in relns2
            if pair1 == pair2 @goto next_pair1 end
        end
        return false
        @label next_pair1
    end
    return true
end

"""
permute_multiplication_table(permutation::Vector{Int64}, multable::Matrix{Int64}) -> Matrix{Int64}
rearrange_multiplication_table(multable::Matrix{Int64}) -> Matrix{Int64}
"""
function permute_multiplication_table(permutation::Vector{Int64}, multable::Matrix{Int64})::Matrix{Int64}
    new_multable = [permutation[ele] for ele in multable]
    new_multable = rearrange_multiplication_table(new_multable)
    return new_multable
end
function rearrange_multiplication_table(multable::Matrix{Int64})::Matrix{Int64}
   result = []
    for i = 1:size(multable, 1)
        for j = 1:size(multable, 1)
            if i == multable[1, j]
                append!(result, multable[:,j ])
                break
            end
        end
    end
    multable = reshape(result, size(multable))
    result = []
    for i = 1:size(multable, 1)
        for j = 1:size(multable, 1)
            if i == multable[j, 1]
                append!(result, multable[j, :])
                break
            end
        end
    end
    multable = transpose(reshape(result, size(multable)))
    return multable
end

end #module GroupWithMultiplicationTable
