include("../math-tool/mathtool.jl")
include("../representation-tool/kspace.jl")
include("../symmetry-operation-tool/pointsymmetry.jl")
include("../symmetry-operation-tool/spinpointsymmetry.jl")
include("../symmetry-group-tool/pointgroup.jl")
include("../representation-tool/grouprepresentation.jl")
include("../representation-tool/groupwithmultiplicationtable.jl")

######################################################################################################################
##################################################FullSpinPointGroup##################################################
######################################################################################################################
module SpinPointGroup

export spg_encode_symmetry
export spg_decode_symmetry
export colAFM_encode_symmetry
export colAFM_decode_symmetry

include("../database/spg_database.jl")
##################################################SpinPointSymmetry###################################################
######################################################################################################################
module SpinPointSymmetry
using LinearAlgebra
using ...MathTool
using ...PointSymmetry
using ..SpinPointGroupDataBase
import ...SpinPointSymmetry as sptsymmetry

export change_symmetry_from_spin_basis_to_magnon_basis
export get_inverse_of_symmetry
export get_symmetry_order
export get_notation_of_spin_point_symmetry
export get_isometry_type
export get_isometry_axis
export get_sense_of_rotation
export is_symmetry_unitary
export is_symmetry_same
export composite_symmetries
export conjugate_symmetry

"""
change_symmetry_from_spin_basis_to_magnon_basis(sym::spg_symtype) -> spg_magnonbasis_symtype
only suitable for magnon in collinear AFM

the final matrix is in basis of local reference frame, where z ∥ collinear spin arrangement
INPUT

Formula
    global reference frame
        Jx, Jy, Jz
    local reference frame (Sz ∥ spin)
        Sx, Sy, Sz
    operator for spin-up and spin-down
        a^+_↑ =  Sx+iSy
        a^+_↓ = -Sx+iSy
    R: rotation matrix between global global reference frame to local reference frame
        (Jx, Jy, Jz) = (Sx, Sy, Sz) R
        if spin ∥ Jx, R = [ 0  0  1; 0  1  0;-1  0  0]
        if spin ∥ Jy, R = [ 0 -1  0; 0  0  1;-1  0  0]
        if spin ∥ Jz, R = [ 1  0  0; 0  1  0; 0  0  1]

    s_global: matrix in spin space of global reference frame
    s_local: matrix in spin space of global reference frame
    s_local = R * s_global * R^-1
    s = s_local[1:2,1:2]

    s Sx = s11*Sx+s12*Sy
      Sy   s21*Sx+s22*Sy

    unitary
    s a^+_↑ = 0.5*( s11-is12+is21+s22)*a^+_↑ + 0.5*(-s11-is12-is21+s22)*a^+_↓
      a^+_↓ = 0.5*(-s11+is12+is21+s22)*a^+_↑ + 0.5*( s11+is12-is21+s22)*a^+_↓
    antiunitary
    s a^+_↑ = 0.5*( s11-is12-is21-s22)*a^+_↑ + 0.5*(-s11-is12+is21-s22)*a^+_↓
      a^+_↓ = 0.5*(-s11+is12-is21-s22)*a^+_↑ + 0.5*( s11+is12+is21-s22)*a^+_↓
Procedure
1. change spin symmetry from global reference frame (s_global) to local reference frame (s_local)
"""
function change_symmetry_from_spin_basis_to_magnon_basis(sym::spg_symtype; SO2::Int64)::spg_magnonbasis_symtype
    if SO2 == 1
        dir = [1,0,0]
    elseif SO2 == 2
        dir = [0,1,0]
    elseif SO2 == 3
        dir = [0,0,1]
    end
    nsym = sptsymmetry.change_symmetry_from_spin_basis_to_magnon_basis((r=sym.r, s=Float64.(sym.s)), SO2=dir)
    nsym = (r=nsym.r, s=Int64.(nsym.s), T=nsym.T)
    return nsym

    rot = [
        [ 0  0  1; 0  1  0;-1  0  0],
        [ 0 -1  0; 0  0  1;-1  0  0],
        [ 1  0  0; 0  1  0; 0  0  1]
    ]
    s_global = deepcopy(sym.s)
    s_local = rot[SO2]*s_global*inv(rot[SO2])
    s = s_local[1:2,1:2]
    T = Int64(round.(det(s_local), digits=5))
    if T == 1
        ns = [
            0.5*( s[1,1]-im*s[1,2]+im*s[2,1]+s[2,2]) 0.5*(-s[1,1]-im*s[1,2]-im*s[2,1]+s[2,2])
            0.5*(-s[1,1]+im*s[1,2]+im*s[2,1]+s[2,2]) 0.5*( s[1,1]+im*s[1,2]-im*s[2,1]+s[2,2])
        ]
    else
        ns = [
            0.5*( s[1,1]-im*s[1,2]-im*s[2,1]-s[2,2]) 0.5*(-s[1,1]-im*s[1,2]+im*s[2,1]-s[2,2])
            0.5*(-s[1,1]+im*s[1,2]-im*s[2,1]-s[2,2]) 0.5*( s[1,1]+im*s[1,2]+im*s[2,1]-s[2,2])
        ]
    end
    nsym = (r=sym.r, s=Int64.(ns), T=T)
    return nsym
end

"""
get_inverse_of_symmetry(encoded::spg_encodetype; encode::Bool)
get_inverse_of_symmetry(sym::spg_symtype; encode::Bool)
"""
function get_inverse_of_symmetry(encoded::Union{spg_encodetype, spg_magnonbasis_encodetype}; encode::Bool)
    sym = spg_decode_symmetry(encoded)
    return get_inverse_of_symmetry(sym, encode=encode)
end
function get_inverse_of_symmetry(sym::spg_symtype; encode::Bool)
    sym_inv = (r=get_inverse_of_matrix(sym.r), s=get_inverse_of_matrix(sym.s))
    if encode
        return spg_encode_symmetry(sym_inv)
    else
        return sym_inv
    end
end
function get_inverse_of_symmetry(sym::spg_magnonbasis_symtype; encode::Bool)
    sym_inv = (r=get_inverse_of_matrix(sym.r), s=get_inverse_of_matrix(sym.s), T=sym.T)
    if encode
        return spg_encode_symmetry(sym_inv)
    else
        return sym_inv
    end
end
#=
function get_inverse_of_symmetry(encoded::colAFM_encodetype[1]; encode::Bool)
    sym = colAFM_decode_symmetry(encoded)
    return get_inverse_of_symmetry(sym, encode=encode)
end
function get_inverse_of_symmetry(sym::colAFM_symtype[1]; encode::Bool)
    sym_inv = (r=get_inverse_of_matrix(sym.r), s=get_inverse_of_matrix(sym.s), T=sym.T)
    if encode
        return colAFM_encode_symmetry(sym_inv)
    else
        return sym_inv
    end
end
=#

"""
get_symmetry_order(encoded::spg_encodetype) -> Int64
get_symmetry_order(sym::spg_symtype) -> Int64
#get order n of a symmetry operation, sym^n=I
"""
function get_symmetry_order(encoded::spg_encodetype)::Int64
    sym = spg_decode_symmetry(encoded)
    return get_symmetry_order(sym)
end
function get_symmetry_order(sym::spg_symtype)::Int64
    order1 = get_symmetry_order_of_point_symmetry(sym.r)
    order2 = get_symmetry_order_of_point_symmetry(sym.s)
    order = lcm(order1, order2)
    return order
end

"""
get_notation_of_spin_point_symmetry(encoded::spg_encodetype)
get_notation_of_spin_point_symmetry(sym::spg_symtype)
"""
function get_notation_of_spin_point_symmetry(encoded::spg_encodetype)
    sym = spg_decode_symmetry(encoded)
    r = get_notation_of_point_symmetry(sym.r)
    s = get_notation_of_point_symmetry(sym.s)
    notation = string("{", s, "||", r, "}")
    return notation
end
function get_notation_of_spin_point_symmetry(sym::spg_symtype)
    r = get_notation_of_point_symmetry(sym.r)
    s = get_notation_of_point_symmetry(sym.s)
    notation = string("{", s, "||", r, "}")
    return notation
end
function get_notation_of_spin_point_symmetry(encoded::spg_magnonbasis_encodetype)
    sym = spg_decode_symmetry(encoded)
    r = get_notation_of_point_symmetry(sym.r)
    if (encoded[2], encoded[3]) == (68, 1)
        s = "E"
    elseif (encoded[2], encoded[3]) == (68, -1)
        s = "dCxT"#"dC_[1, 0, 0](π)T"
    elseif (encoded[2], encoded[3]) == (66, 1)
        s = "Σz"
    elseif (encoded[2], encoded[3]) == (66, -1)
        s = "dΣzCxT"#"dΣzC_[1, 0, 0](π)T"
    elseif (encoded[2], encoded[3]) == (12, 1)
        s = "dE"
    elseif (encoded[2], encoded[3]) == (12, -1)
        s = "CxT"#"C_[1, 0, 0](π)T"
    elseif (encoded[2], encoded[3]) == (14,  1)
        s = "dΣz"
    elseif (encoded[2], encoded[3]) == (14, -1)
        s = "ΣzCxT"#"ΣzC_[1, 0, 0](π)T"
    elseif (encoded[2], encoded[3]) == (52,  1)
        s = "dCx"#"dC_[1, 0, 0](π)"
    elseif (encoded[2], encoded[3]) == (52, -1)
        s = "T"
    elseif (encoded[2], encoded[3]) == (34, 1)
        s = "dΣzCx"
    elseif (encoded[2], encoded[3]) == (34, -1)
        s = "ΣzT"
    elseif (encoded[2], encoded[3]) == (28,  1)
        s = "Cx"
    elseif (encoded[2], encoded[3]) == (28, -1)
        s = "dT"
    elseif (encoded[2], encoded[3]) == (46,  1)
        s = "ΣzCx"
    elseif (encoded[2], encoded[3]) == (46, -1)
        s = "dΣzT"
    else
        println((encoded[2], encoded[3]))
    end
    #s = get_notation_of_point_symmetry(sym.s)
    notation = string("{", s, "||", r, "}")
    return notation
end
function get_notation_of_spin_point_symmetry(sym::spg_magnonbasis_symtype)
    encoded = spg_encode_symmetry(sym)
    return get_notation_of_spin_point_symmetry(encoded)
end

"""
get_isometry_type(encoded::spg_encodetype) -> Tuple{Int64, Int64}
get_isometry_type(sym::spg_symtype) -> Tuple{Int64, Int64}
isometry_type
operation    -6 -4 -3 -2 -1  1  2  3  4  6
trace        -2 -1  0  1 -3  3 -1  0  1  2
determinant  -1 -1 -1 -1 -1  1  1  1  1  1
index         1  2  3  4  5  6  7  8  9 10
"""
function get_isometry_type(encoded::spg_encodetype)::Tuple{Int64, Int64}
    sym = spg_decode_symmetry(encoded)
    return get_isometry_type(sym)
end
function get_isometry_type(sym::spg_symtype)::Tuple{Int64, Int64}
    r = get_isometry_type_of_point_symmetry(sym[1])
    s = get_isometry_type_of_point_symmetry(sym[2])
    if r == nothing
        println(1,sym)
    elseif s == nothing
        println(2,sym)
    end
    iso_type = (r, s)
    return iso_type
end
#=
function get_isometry_type(encoded::colAFM_encodetype[1])::Tuple{Int64, Int64}
    sym = colAFM_decode_symmetry(encoded)
    return get_isometry_type(sym)
end
function get_isometry_type(sym::colAFM_symtype[1])::Tuple{Int64, Int64}
    r = get_isometry_type_of_point_symmetry(sym.r)
    encoded = colAFM_encode_symmetry(sym)
    if encoded[2] in [68, 66, 12, 14]
        if sym.T == 1
            s = 6
        else
            s = 4
        end
    elseif encoded[2] in [52, 34, 28, 46]
        if sym.T == 1
            s = 7
        else
            s = 5
        end
    else
        error("error occurs at get_isometry_type")
    end
    iso_type = (r, s)
    return iso_type
end
=#

"""
get_isometry_axis(encoded::spg_encodetype) -> Tuple{Vector{Int64}, Vector{Int64}}
get_isometry_axis(sym::spg_symtype) -> Tuple{Vector{Int64}, Vector{Int64}}
#get isometry axis of a symmetry operation
"""
function get_isometry_axis(encoded::spg_encodetype)::Tuple{Vector{Int64}, Vector{Int64}}
    sym = spg_decode_symmetry(encoded)
    iso_axis = (get_isometry_axis_of_point_symmetry(sym.r), get_isometry_axis_of_point_symmetry(sym.s))
    return iso_axis
end
function get_isometry_axis(sym::spg_symtype)::Tuple{Vector{Int64}, Vector{Int64}}
    iso_axis = (get_isometry_axis_of_point_symmetry(sym.r), get_isometry_axis_of_point_symmetry(sym.s))
    return iso_axis
end

"""
get_sense_of_rotation(encoded::spg_encodetype) -> Int64
get_sense_of_rotation(sym::spg_symtype) -> Int64
#get sense of rotation of a symmetry operation
"""
function get_sense_of_rotation(encoded::spg_encodetype)::Int64
    sym = spg_decode_symmetry(encoded)
    rot_sense = (get_sense_of_rotation_of_point_symmetry(sym.r), get_sense_of_rotation_of_point_symmetry(sym.s))
    return rot_sense
end
function get_sense_of_rotation(sym::spg_symtype)::Int64
    rot_sense = (get_sense_of_rotation_of_point_symmetry(sym.r), get_sense_of_rotation_of_point_symmetry(sym.s))
    return rot_sense
end

"""
is_symmetry_unitary(encoded::spg_encodetype) -> Bool
is_symmetry_unitary(sym::spg_symtype) -> Bool
"""
function is_symmetry_unitary(encoded::spg_encodetype)::Bool
    sym = spg_decode_symmetry(encoded)
    return is_symmetry_unitary(sym)
end
function is_symmetry_unitary(sym::spg_symtype)::Bool
    s = sym.s
    if Int64(det(s)) == 1
        return true
    elseif Int64(det(s)) == -1
        return false
    else
        error("Error occurs at is_symmetry_unitary: get wrong determinant of spin part operation")
    end
end
function is_symmetry_unitary(encoded::spg_magnonbasis_encodetype)::Bool
    if encoded[3] == 1
        return true
    else
        return false
    end
end
function is_symmetry_unitary(sym::spg_magnonbasis_symtype)::Bool
    if sym.T == 1
        return true
    else
        return false
    end
end

"""
is_symmetry_same(sym1::Union{spg_symtype, spg_magnonbasis_symtype}, encoded2::Union{spg_symtype, spg_magnonbasis_symtype}) -> Bool
is_symmetry_same(encoded1::Union{spg_encodetype, spg_magnonbasis_encodetype}, encoded2::Union{spg_encodetype, spg_magnonbasis_encodetype}) -> Bool
"""
function is_symmetry_same(sym1::Union{spg_symtype, spg_magnonbasis_symtype}, sym2::Union{spg_symtype, spg_magnonbasis_symtype})::Bool
    encoded1 = spg_encode_symmetry(sym1)
    encoded2 = spg_encode_symmetry(sym2)
    return is_symmetry_same(encoded1, encoded2)
end
function is_symmetry_same(encoded1::Union{spg_encodetype, spg_magnonbasis_encodetype}, encoded2::Union{spg_encodetype, spg_magnonbasis_encodetype})::Bool
    if encoded1 == encoded2
        return true
    else
        return false
    end
end

"""
composite_symmetries(encoded::Union{spg_encodetype, spg_magnonbasis_encodetype}; order::Integer, encode::Bool)::Union{spg_encodetype, spg_symtype, spg_magnonbasis_symtype, spg_magnonbasis_encodetype}
composite_symmetries(sym::Union{spg_symtype, spg_magnonbasis_symtype}; order::Integer, encode::Bool)::Union{spg_encodetype, spg_symtype, spg_magnonbasis_symtype, spg_magnonbasis_encodetype}
composite_symmetries(encodeds::Union{spg_encodetype, spg_magnonbasis_encodetype}...; encode::Bool)::Union{spg_encodetype, spg_symtype, spg_magnonbasis_symtype, spg_magnonbasis_encodetype}
composite_symmetries(syms::Union{spg_symtype, spg_magnonbasis_symtype}...; encode::Bool)::Union{spg_encodetype, spg_symtype, spg_magnonbasis_symtype, spg_magnonbasis_encodetype}
"""
function composite_symmetries(encoded::Union{spg_encodetype, spg_magnonbasis_encodetype}; order::Integer, encode::Bool)::Union{spg_encodetype, spg_symtype, spg_magnonbasis_symtype, spg_magnonbasis_encodetype}
    sym = spg_decode_symmetry(encoded)
    return composite_symmetries(sym; order=order, encode=encode)
end
function composite_symmetries(sym::Union{spg_symtype, spg_magnonbasis_symtype}; order::Integer, encode::Bool)::Union{spg_encodetype, spg_symtype, spg_magnonbasis_symtype, spg_magnonbasis_encodetype}
    if order == 1
        if encode
            return spg_encode_symmetry(sym)
        else
            return sym
        end
    elseif order > 1
        return composite_symmetries([sym for i = 1:order]..., encode=encode)
    else
        error("Error occurs at composite_symmetries: order should be a positive integer")
    end
end
function composite_symmetries(encodeds::Union{spg_encodetype, spg_magnonbasis_encodetype}...; encode::Bool)::Union{spg_encodetype, spg_symtype, spg_magnonbasis_symtype, spg_magnonbasis_encodetype}
    syms = (spg_decode_symmetry(encoded) for encoded in encodeds)
    return composite_symmetries(syms..., encode=encode)
end
function composite_symmetries(syms::Union{spg_symtype, spg_magnonbasis_symtype}...; encode::Bool)::Union{spg_encodetype, spg_symtype, spg_magnonbasis_symtype, spg_magnonbasis_encodetype}
    if length(syms) == 1
        if encode
            return spg_encode_symmetry(syms[1])
        else
            return syms[1]
        end
    elseif length(syms) == 2
        if isa(syms[1], spg_symtype)
            nsym = (r=syms[1].r*syms[2].r, s=syms[1].s*syms[2].s)
        else
            nsym = (r=syms[1].r*syms[2].r, s=syms[1].s*syms[2].s, T=syms[1].T*syms[2].T)
        end
        if encode
            return spg_encode_symmetry(nsym)
        else
            return nsym
        end
    else
        return composite_symmetries(composite_symmetries(syms[1], syms[2], encode=false), syms[3:end]..., encode=encode)
    end
end

"""
conjugate_symmetry(sym_encoded::spg_encodetype, conjg_sym_encoded::spg_encodetype; encode::Bool) -> Union{spg_encodetype, spg_symtype}
conjugate_symmetry(sym::spg_symtype, conjg_sym::spg_symtype; encode::Bool) -> Union{spg_encodetype, spg_symtype}
conjugate_symmetry(sym_encoded::colAFM_encodetype[1], conjg_sym_encoded::colAFM_encodetype[1]; encode::Bool) -> Union{colAFM_encodetype[1], colAFM_symtype[1]}
conjugate_symmetry(sym::colAFM_symtype[1], conjg_sym::colAFM_symtype[1]; encode::Bool) -> Union{colAFM_encodetype[1], colAFM_symtype[1]}
#matrices are in the same coordinate systems by default
#conjg_sym^{-1}*sym*conjg_sym
"""
function conjugate_symmetry(sym_encoded::Union{spg_encodetype, spg_magnonbasis_encodetype}, conjg_sym_encoded::Union{spg_encodetype, spg_magnonbasis_encodetype}; encode::Bool)::Union{spg_encodetype, spg_symtype, spg_magnonbasis_symtype, spg_magnonbasis_encodetype}
    sym = spg_decode_symmetry(sym_encoded)
    conjg_sym = spg_decode_symmetry(conjg_sym_encoded)
    return conjugate_symmetry(sym, conjg_sym, encode=encode)
end
function conjugate_symmetry(sym::spg_symtype, conjg_sym::spg_symtype; encode::Bool)::Union{spg_encodetype, spg_symtype}
    return composite_symmetries((r=Int64.(inv(conjg_sym.r)), s=Int64.(inv(conjg_sym.s))), sym, conjg_sym, encode=encode)
end
function conjugate_symmetry(sym::spg_magnonbasis_symtype, conjg_sym::spg_magnonbasis_symtype; encode::Bool)::Union{spg_magnonbasis_encodetype, spg_magnonbasis_symtype}
    return composite_symmetries((r=Int64.(inv(conjg_sym.r)), s=Int64.(inv(conjg_sym.s)), T=conjg_sym.T), sym, conjg_sym, encode=encode)
end

end #(moduleSpinPointSymmetry)

###################################################GroupInformation###################################################
######################################################################################################################
module GroupInformation
using LinearAlgebra
using ..SpinPointGroupDataBase
using ..SpinPointSymmetry
import ...PointGroup as pg
import ...SpinPointSymmetry as sptsymmetry
import ...GroupWithMultiplicationTable as gwmt

export get_multiplication_table_of_nontrivial_spin_point_group
export get_multiplication_table_of_spin_point_group
export get_multiplication_table_of_group
export isa_magnetic_point_group
export get_full_conjugacy_classes_of_nontrivial_spin_point_group
export get_full_conjugacy_classes_of_spin_point_group
export get_full_conjugacy_classes_of_group
export is_group_same
export isa_subgroup
export get_full_subgroups_of_nontrivial_spin_point_group
export get_full_subgroups_of_spin_point_group
export get_full_subgroups_of_group
export is_nontrivial_spin_point_group_unitary
export is_spin_point_group_unitary
export is_group_unitary
export get_maximal_unitary_subgroup_of_nontrivial_spin_point_group
export get_maximal_unitary_subgroup_of_spin_point_group
export get_maximal_unitary_subgroup
export get_maximal_antiunitary_set_of_nontrivial_spin_point_group
export get_maximal_antiunitary_set_of_spin_point_group
export get_maximal_antiunitary_set
export get_class_of_nontrivial_spin_point_group_of_full_spin_point_group
export get_class_of_nontrivial_spin_point_group
export get_class_of_spin_point_group
export get_notation_of_nontrivial_spin_point_group
export get_index_of_nontrivial_spin_point_group_by_notation
export get_isometry_type_table
export get_full_symmetry_operations_of_nontrivial_spin_point_group
export get_full_symmetry_operations_of_spin_point_group
export get_reduced_group
export isa_set
export isa_group
export generate_full_spin_point_group_of_nontrivial_spin_point_group

"""
get_multiplication_table_of_nontrivial_spin_point_group(nspgidx::Int64) -> Matrix{Int64}
"""
function get_multiplication_table_of_nontrivial_spin_point_group(nspgidx::Int64)::Matrix{Int64}
    syms = get_full_symmetry_operations_of_nontrivial_spin_point_group(nspgidx, encode=false)
    return get_multiplication_table_of_nontrivial_spin_point_group(syms)
end

"""
get_multiplication_table_of_spin_point_group(spgidx::Int64) -> Matrix{Int64}
"""
function get_multiplication_table_of_spin_point_group(spgidx::Int64)::Matrix{Int64}
    syms = get_full_symmetry_operations_of_spin_point_group(spgidx, encode=false)
    return get_multiplication_table(syms)
end

"""
get_multiplication_table_of_group(encodeds::Vector{spg_encodetype}) -> Matrix{Int64}
get_multiplication_table_of_group(syms::Vector{spg_symtype, spg_magnonbasis_symtype})::Matrix{Int64}
"""
function get_multiplication_table_of_group(encodeds::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}})::Matrix{Int64}
    syms = [spg_decode_symmetry(encoded) for encoded in encodeds]
    return get_multiplication_table_of_group(syms)
end
function get_multiplication_table_of_group(syms::Vector{<:Union{spg_symtype, spg_magnonbasis_symtype}})::Matrix{Int64}
    if !isa_group(syms) error("Error occurs at get_multiplication_table_of_group: symmetry operations don't form a group") end
    order = length(syms)
    mul_table = zeros(Int64, order, order)
    for i = 1:order
        for j = 1:order
            nsym = composite_symmetries(syms[i], syms[j], encode=false)
            idx = findall(sym->is_symmetry_same(sym, nsym), syms)[1]
            mul_table[i, j] = idx
        end
    end
    return mul_table
end

"""
isa_magnetic_point_group(nspgnum::Int64) -> Bool
isa_magnetic_point_group(syms::Vector{spg_symtype}; tor=1e-5) -> Bool
"""
function isa_magnetic_point_group(nspgnum::Int64)::Bool
    syms = get_full_symmetry_operations_of_nontrivial_spin_point_group(nspgnum, encode=false)
    return isa_magnetic_point_group(syms)
end
function isa_magnetic_point_group(syms::Vector{spg_symtype}; tor=1e-5)::Bool
    for sym in syms
        r = sym.r * det(sym.r)
        s = sym.s * det(sym.s)
        if !all(abs.(round.(r-s, digits=5)) .< tor) return false end
    end
    return true
end

"""
get_full_conjugacy_classes_of_nontrivial_spin_point_group(nspgidx::Int64; encode::Bool) -> Union{Vector{Vector{spg_encodetype}}, Vector{Vector{spg_symtype}}}
"""
function get_full_conjugacy_classes_of_nontrivial_spin_point_group(nspgidx::Int64; encode::Bool)::Union{Vector{Vector{spg_encodetype}}, Vector{Vector{spg_symtype}}}
    syms = get_full_symmetry_operations_of_nontrivial_spin_point_group(nspgidx, encode=true)
    return get_full_conjugacy_classes_of_group(syms, encode=encode)
end

"""
get_full_conjugacy_classes_of_spin_point_group(spgidx::Int64; encode::Bool) -> Union{Vector{Vector{spg_encodetype}}, Vector{Vector{spg_symtype}}}
"""
function get_full_conjugacy_classes_of_spin_point_group(spgidx::Int64; encode::Bool)::Union{Vector{Vector{spg_encodetype}}, Vector{Vector{spg_symtype}}}
    syms = get_full_symmetry_operations_of_spin_point_group(spgidx, encode=true)
    return get_full_conjugacy_classes_of_group(syms, encode=encode)
end

"""
get_full_conjugacy_classes_of_group(syms::Vector{spg_symtype}; encode::Bool) -> Union{Vector{Vector{spg_encodetype}}, Vector{Vector{spg_symtype}}}
get_full_conjugacy_classes_of_group(syms::Vector{spg_encodetype}; encode::Bool) -> Union{Vector{Vector{spg_encodetype}}, Vector{Vector{spg_symtype}}}
"""
function get_full_conjugacy_classes_of_group(syms::Vector{spg_symtype}; encode::Bool)::Union{Vector{Vector{spg_encodetype}}, Vector{Vector{spg_symtype}}}
    encodeds = [spg_encode_symmetry(sym) for sym in syms]
    return get_full_conjugacy_classes_of_group(encodeds, encode=encode)
end
function get_full_conjugacy_classes_of_group(syms::Vector{spg_encodetype}; encode::Bool)::Union{Vector{Vector{spg_encodetype}}, Vector{Vector{spg_symtype}}}
    if !isa_group(syms) error("Error occurs at get_full_conjugacy_classes_of_group: symmetry operations don't form a group") end
    syms_remain = [sym for sym in syms if sym !== (16484, 16484)]
    conjg_classes = Vector{Vector{spg_encodetype}}([[(16484, 16484)]])
    while syms_remain != []
        class = Vector{spg_encodetype}()
        push!(class, syms_remain[1])
        for sym in syms
            nsym = conjugate_symmetry(syms_remain[1], sym, encode=true)
            if !(nsym in syms) error("Error occurs at get_conjugacy_classes_of_point_groups: group not closed") end
            if !(nsym in class) push!(class, nsym) end
        end
        push!(conjg_classes, class)
        syms_remain = [sym for sym in syms_remain if !(sym in class)]
    end
    return conjg_classes
end

"""
is_group_same(group1::Vector{spg_symtype}, group2::Vector{spg_encodetype}) -> Bool
is_group_same(group1::Vector{spg_encodetype}, group2::Vector{spg_symtype}) -> Bool
is_group_same(group1::Vector{spg_encodetype}, group2::Vector{spg_encodetype}) -> Bool
#check if two groups are exactly the same without any isomorphism
"""
function is_group_same(group1::Vector{spg_symtype}, group2::Vector{spg_encodetype})::Bool
    group1 = [spg_encode_symmetry(sym) for sym in group1]
    return is_group_same(group1, group2)
end
function is_group_same(group1::Vector{spg_encodetype}, group2::Vector{spg_symtype})::Bool
    group2 = [spg_encode_symmetry(sym) for sym in group2]
    return is_group_same(group1, group2)
end
function is_group_same(group1::Vector{spg_encodetype}, group2::Vector{spg_encodetype})::Bool
    if length(group1) != length(group2) return false end
    for sym in group1
        if !(sym in group2) return false end
    end
    return true
end

"""
#isa_subgroup(subnum::Int64, groupnum::Int64) -> Bool
isa_subgroup(group1::Vector{spg_encodetype}, group2::Vector{spg_encodetype}) -> Bool
"""
#=
function isa_subgroup(subnum::Int64, groupnum::Int64)::Bool
    group1 = get_full_symmetry_operations_of_nontrivial_spin_point_group(subnum, encode=true)
    group2 = get_full_symmetry_operations_of_nontrivial_spin_point_group(groupnum, encode=true)
    return isa_subgroup(group1, group2)
end
=#
function isa_subgroup(group1::Vector{spg_encodetype}, group2::Vector{spg_encodetype})::Bool
    for encoded in group1
        if !(encoded in group2) return false end
    end
    return true
end

"""
get_full_subgroups_of_nontrivial_spin_point_group(nspgnum::Int64; encode::Bool) -> Union{Vector{Vector{spg_encodetype}}, Vector{Vector{spg_symtype}}}
#get all full subgroups of given group from NontrivialSpinPointGroupDataBase
"""
function get_full_subgroups_of_nontrivial_spin_point_group(nspgidx::Int64; encode::Bool)::Union{Vector{Vector{spg_encodetype}}, Vector{Vector{spg_symtype}}}
    if nspgidx < 1 || nspgidx > 598 error("Error occurs at get_full_subgroups_of_nontrivial_spin_point_group: nontrivial spion point group number should be in range of [1, 598]") end
    subgroups = nspg_full_subgroups[nspgidx]
    if encode
        subgroups = [collect(subgroup) for subgroup in subgroups]
    else
        subgroups = [[spg_decode_symmetry(sym) for sym in subgroup] for subgroup in subgroups]
    end
    return subgroups
end


"""
get_full_subgroups_of_spin_point_group(spgidx::Int64; encode::Bool) -> Union{Vector{Vector{spg_encodetype}}, Vector{Vector{spg_symtype}}}
#get all full subgroups of given group from NontrivialSpinPointGroupDataBase
"""
function get_full_subgroups_of_spin_point_group(spgidx::Int64; encode::Bool)::Union{Vector{Vector{spg_encodetype}}, Vector{Vector{spg_symtype}}}
    if spgidx < 1 || spgidx > 1263 error("Error occurs at get_full_subgroups_of_spin_point_group: spin point group number should be in range of [1, 941]") end
    subgroups = spg_full_subgroups[spgidx]
    if encode
        subgroups = [collect(subgroup) for subgroup in subgroups]
    else
        subgroups = [[spg_decode_symmetry(sym) for sym in subgroup] for subgroup in subgroups]
    end
    return subgroups
end

"""
get_full_subgroups_of_group(encodeds::Vector{spg_encodetype}; encode::Bool) -> Union{Vector{Vector{spg_encodetype}}, Vector{Vector{spg_symtype}}}
"""
function get_full_subgroups_of_group(encodeds::Vector{spg_encodetype}; encode::Bool)::Union{Vector{Vector{spg_encodetype}}, Vector{Vector{spg_symtype}}}
    #=get all unitary subgroups=#
    uencodeds = get_maximal_unitary_subgroup(encodeds, encode=true)
    nspgidx = get_class_of_nontrivial_spin_point_group(uencodeds)
    nspgencodeds = get_full_symmetry_operations_of_nontrivial_spin_point_group(nspgidx, encode=true)
    nspgsubs = get_full_subgroups_of_nontrivial_spin_point_group(nspgidx, encode=true)
    if is_group_same(uencodeds, nspgencodeds)
        usubs = nspgsubs
    else
        multable1 = get_multiplication_table_of_group(nspgencodeds)
        multable2 = get_multiplication_table_of_group(uencodeds)
        perm = gwmt.get_isomorphism(multable1, multable2)
        usubs = [[uencodeds[perm[findall(x->x==encoded, nspgencodeds)[1]]] for encoded in sub] for sub in nspgsubs]
    end

    #=get all antiunitary subgroups=#
    aencodeds = get_maximal_antiunitary_set(encodeds, encode=true)
    asubs = []
    for sub in usubs
        for a in aencodeds
            group = [sub; [composite_symmetries(a, encoded, encode=true) for encoded in sub]]
            if !isa_group(group) continue end
            for sub in asubs
                if is_group_same(group, sub) @goto next end
            end
            push!(asubs, group)
            @label next
        end
    end

    subs = Vector{spg_encodetype}[usubs; asubs]
    if !encode subs = Vector{spg_symtype}[[spg_decode_symmetry(encoded) for encoded in sub] for sub in subs] end
    return subs
end

"""
is_nontrivial_spin_point_group_unitary(nspgidx::Int64) -> Bool
"""
function is_nontrivial_spin_point_group_unitary(nspgidx::Int64)::Bool
    syms = get_full_symmetry_operations_of_nontrivial_spin_point_group(nspgidx, encode=false)
    return is_group_unitary(syms)
end

"""
is_spin_point_group_unitary(spgidx::Int64) -> Bool
"""
function is_spin_point_group_unitary(spgidx::Int64)::Bool
    syms = get_full_symmetry_operations_of_spin_point_group(spgidx, encode=false)
    return is_group_unitary(syms)
end

"""
is_group_unitary(syms::Union{Vector{spg_encodetype}, Vector{spg_symtype}}) -> Bool
"""
function is_group_unitary(syms::Union{Vector{spg_encodetype}, Vector{spg_symtype}})::Bool
    for sym in syms
        if !is_symmetry_unitary(sym) return false end
    end
    return true
end

"""
get_maximal_unitary_subgroup_of_nontrivial_spin_point_group(nspgidx::Int64; encode::Bool) -> Union{Vector{spg_encodetype}, Vector{spg_symtype}}
"""
function get_maximal_unitary_subgroup_of_nontrivial_spin_point_group(nspgidx::Int64; encode::Bool)::Union{Vector{spg_encodetype}, Vector{spg_symtype}}
    syms = get_full_symmetry_operations_of_nontrivial_spin_point_group(nspgidx, encode=false)
    return get_maximal_unitary_subgroup(syms, encode=encode)
end

"""
get_maximal_unitary_subgroup_of_spin_point_group(spgidx::Int64; encode::Bool) -> Union{Vector{spg_encodetype}, Vector{spg_symtype}}
"""
function get_maximal_unitary_subgroup_of_spin_point_group(spgidx::Int64; encode::Bool)::Union{Vector{spg_encodetype}, Vector{spg_symtype}}
    syms = get_full_symmetry_operations_of_spin_point_group(spgidx, encode=false)
    return get_maximal_unitary_subgroup(syms, encode=encode)
end

"""
get_maximal_unitary_subgroup(syms::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}}; encode::Bool) -> Vector{<:Union{spg_encodetype, spg_symtype, spg_magnonbasis_encodetype, spg_magnonbasis_symtype}}
get_maximal_unitary_subgroup(syms::Vector{<:Union{spg_symtype, spg_magnonbasis_symtype}}; encode::Bool) -> Vector{<:Union{spg_encodetype, spg_symtype, spg_magnonbasis_encodetype, spg_magnonbasis_symtype}}
"""
function get_maximal_unitary_subgroup(syms::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}}; encode::Bool)::Vector{<:Union{spg_encodetype, spg_symtype, spg_magnonbasis_encodetype, spg_magnonbasis_symtype}}
    unitary_subgroup = [sym for sym in syms if is_symmetry_unitary(sym)]
    if !encode unitary_subgroup = [spg_decode_symmetry(sym) for sym in unitary_subgroup] end
    return unitary_subgroup
end
function get_maximal_unitary_subgroup(syms::Vector{<:Union{spg_symtype, spg_magnonbasis_symtype}}; encode::Bool)::Vector{<:Union{spg_encodetype, spg_symtype, spg_magnonbasis_encodetype, spg_magnonbasis_symtype}}
    unitary_subgroup = [sym for sym in syms if is_symmetry_unitary(sym)]
    if encode unitary_subgroup = [spg_encode_symmetry(sym) for sym in unitary_subgroup] end
    return unitary_subgroup
end

"""
get_maximal_antiunitary_set_of_nontrivial_spin_point_group(nspgidx::Int64; encode::Bool) -> Union{Vector{spg_encodetype}, Vector{spg_symtype}}
"""
function get_maximal_antiunitary_set_of_nontrivial_spin_point_group(nspgidx::Int64; encode::Bool)::Union{Vector{spg_encodetype}, Vector{spg_symtype}}
    syms = get_full_symmetry_operations_of_nontrivial_spin_point_group(nspgidx, encode=false)
    return get_maximal_antiunitary_set(syms, encode=encode)
end

"""
get_maximal_antiunitary_set_of_spin_point_group(spgidx::Int64; encode::Bool) -> Union{Vector{spg_encodetype}, Vector{spg_symtype}}
"""
function get_maximal_antiunitary_set_of_spin_point_group(spgidx::Int64; encode::Bool)::Union{Vector{spg_encodetype}, Vector{spg_symtype}}
    syms = get_full_symmetry_operations_of_spin_point_group(nspgidx, encode=false)
    return get_maximal_antiunitary_set(syms, encode=encode)
end

"""
get_maximal_antiunitary_set(syms::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}}; encode::Bool) -> Vector{<:Union{spg_encodetype, spg_symtype, spg_magnonbasis_encodetype, spg_magnonbasis_symtype}}
get_maximal_antiunitary_set(syms::Vector{<:Union{spg_symtype, spg_magnonbasis_symtype}}; encode::Bool) -> Vector{<:Union{spg_encodetype, spg_symtype, spg_magnonbasis_encodetype, spg_magnonbasis_symtype}}
"""
function get_maximal_antiunitary_set(syms::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}}; encode::Bool)::Vector{<:Union{spg_encodetype, spg_symtype, spg_magnonbasis_encodetype, spg_magnonbasis_symtype}}
    antiunitary_set = [sym for sym in syms if !is_symmetry_unitary(sym)]
    if !encode antiunitary_set = [spg_decode_symmetry(sym) for sym in antiunitary_set] end
    return antiunitary_set
end
function get_maximal_antiunitary_set(syms::Vector{<:Union{spg_symtype, spg_magnonbasis_symtype}}; encode::Bool)::Vector{<:Union{spg_encodetype, spg_symtype, spg_magnonbasis_encodetype, spg_magnonbasis_symtype}}
    antiunitary_set = [sym for sym in syms if !is_symmetry_unitary(sym)]
    if encode antiunitary_set = [spg_encode_symmetry(sym) for sym in antiunitary_set] end
    return antiunitary_set
end

#=
"""
get_class_of_nontrivial_spin_point_group_of_full_spin_point_group(fspgidx::Int64)::Int64
"""
function get_class_of_nontrivial_spin_point_group_of_full_spin_point_group(fspgidx::Int64)::Int64
    if 1 <= fspgidx <= 941
        nspgidx = nosocfspg_data[fspgidx].nspg
    else
        nspgidx = 0
    end
    return nspgidx
end
=#
"""
get_class_of_nontrivial_spin_point_group(syms::Union{Vector{spg_encodetype}, Vector{spg_symtype}}) -> Int64
"""
function get_class_of_nontrivial_spin_point_group(syms::Union{Vector{spg_encodetype}, Vector{spg_symtype}})::Int64
    iso_type_table = get_isometry_type_table(syms...)
    for idx in eachindex(nspg_data)
        syms_idx = get_full_symmetry_operations_of_nontrivial_spin_point_group(idx, encode=true)
        iso_type_table_idx = get_isometry_type_table(syms_idx...)
        if all(iso_type_table .== iso_type_table_idx) return idx end
    end
    return 0
end

"""
get_class_of_spin_point_group(syms::Union{Vector{spg_encodetype}, Vector{spg_symtype}}; SO2::Int64) -> Int64
"""
function get_class_of_spin_point_group(syms::Union{Vector{spg_encodetype}, Vector{spg_symtype}}; SO2::Int64)::Int64
    if SO2 != 0 syms = get_reduced_operations_of_spin_point_group_with_SO2(syms, SO2=SO2) end
    table = get_isometry_type_table(syms...)
    for idx in eachindex(spg_data)
        if !(SO2 == spg_data[idx].SO2 == 0 || SO2*spg_data[idx].SO2 == 0) || (length(syms) != spg_data[idx].order) continue end
        syms_idx = get_full_symmetry_operations_of_spin_point_group(idx, encode=true)
        table_idx = get_isometry_type_table(syms_idx...)
        if all(table .== table_idx) return idx end
    end
    return 0
end

"""
get_reduced_operations_of_spin_point_group_with_SO2(syms::Vector{spg_symtype}, SO2::Int64) -> Vector{spg_symtype}
"""
function get_reduced_operations_of_spin_point_group_with_SO2(encodeds::Vector{spg_encodetype}; SO2::Int64)::Vector{spg_symtype}
    syms = [spg_decode_symmetry(encoded) for encoded in encodeds]
    return get_reduced_operations_of_spin_point_group_with_SO2(syms, SO2=SO2)
end
function get_reduced_operations_of_spin_point_group_with_SO2(syms::Vector{spg_symtype}; SO2::Int64)::Vector{spg_symtype}
    rsyms = spg_symtype[]
    if SO2 == 1
        so2axis = [1, 0, 0]
    elseif SO2 == 2
        so2axis = [0, 1, 0]
    elseif SO2 == 3
        so2axis = [0, 0, 1]
    end
    for mat in syms
        axis = pg.get_isometry_axis(mat.s)
        if axis == [0, 0, 0]
            push!(rsyms, mat)
        elseif axis == so2axis
            if det(mat.s) > 0
                push!(rsyms, (r=mat.r, s=[1 0 0;0 1 0;0 0 1]))
            else
                push!(rsyms, (r=mat.r, s=[-1 0 0;0 -1 0;0 0 -1]))
            end
        elseif sum(axis.*so2axis) == 0
            if det(mat.s) > 0
                push!(rsyms, (r=mat.r, s=[1 0 0;0 -1 0;0 0 -1]))
            else
                push!(rsyms, (r=mat.r, s=[-1 0 0;0 1 0;0 0 1]))
            end
        else
            error(1)
        end
    end
    return rsyms
end

"""
get_notation_of_nontrivial_spin_point_group(syms::Union{Vector{spg_encodetype}, Vector{spg_symtype}})
get_notation_of_nontrivial_spin_point_group(idx::Int64)
"""
function get_notation_of_nontrivial_spin_point_group(syms::Union{Vector{spg_encodetype}, Vector{spg_symtype}}; subscript=false)
    idx = get_class_of_nontrivial_spin_point_group(syms)
    return get_notation_of_nontrivial_spin_point_group(idx, subscript=subscript)
end
function get_notation_of_nontrivial_spin_point_group(idx::Int64; subscript=false)
    symbol = nspg_data[idx].notation
    if subscript
        i = 1
        while i <= length(symbol)
            if symbol[i] == ')' && symbol[i-1] in ['2', '3', '4', '6', 'm']
                encoded = [encoded for encoded in get_full_symmetry_operations_of_nontrivial_spin_point_group(idx, encode=true) if !(encoded[2] in [16484, 3198])][1]
                axis = get_isometry_axis(encoded)[2]
                if axis == [1, 0, 0]
                    script = "x"
                elseif axis == [0, 1, 0]
                    script = "y"
                elseif axis == [0, 0, 1]
                    script = "z"
                elseif axis == [1, 1, 0]
                    script = "xy"
                elseif axis == [1, 0, 1]
                    script = "xz"
                elseif axis == [0, 1, 1]
                    script = "yz"
                else
                    println(axis)
                    error()
                end
                symbol = string(symbol[1:i-1], script, symbol[i:end])
                i+=2
            else
                i+=1
            end
        end
    end
    return symbol
end

"""
get_index_of_nontrivial_spin_point_group_by_notation(notation) -> Int64
"""
function get_index_of_nontrivial_spin_point_group_by_notation(notation)::Int64
    for idx = 1:598
        if notation == get_notation_of_nontrivial_spin_point_group(idx)
            return idx
        end
    end
    error("Error occurs at get_index_of_nontrivial_spin_point_group_by_notation: can not find index")
end

"""
get_isometry_type_table(syms::Union{Vector{spg_encodetype}, Vector{spg_symtype}}) -> Matrix{Int64}
get_isometry_type_table(encodeds::spg_encodetype...) -> Matrix{Int64}
get_isometry_type_table(syms::spg_symtype...) -> Matrix{Int64}
isometry_type_table
         (spinspace)
         idx  1  2  3  4  5  6  7  8  9 10
         op  -6 -4 -3 -2 -1  1  2  3  4  6
(rspace)
idx op
1   -6
2   -4
3   -3
4   -2
5   -1
6    1
7    2
8    3
9    4
10   6
"""
function get_isometry_type_table(syms::Union{Vector{spg_encodetype}, Vector{spg_symtype}})::Matrix{Int64}
    return get_isometry_type_table(syms...)
end
function get_isometry_type_table(encodeds::spg_encodetype...)::Matrix{Int64}
    syms = [spg_decode_symmetry(encoded) for encoded in encodeds]
    return get_isometry_type_table(syms...)
end
function get_isometry_type_table(syms::spg_symtype...)::Matrix{Int64}
    error = zeros(Int64, 10, 10)
    table = zeros(Int64, 10, 10)
    for sym in syms
        r, s = get_isometry_type(sym)
        if r == 0 || s == 0
            return error
        else
            table[r, s] += 1
        end
    end
    return table
end

"""
get_full_symmetry_operations_of_nontrivial_spin_point_group(nspgidx::Int64; encode::Bool) -> Union{Vector{spg_encodetype}, Vector{spg_symtype}}
given a nspg number, find all its symmetry operations
"""
function get_full_symmetry_operations_of_nontrivial_spin_point_group(nspgidx::Int64; encode::Bool)::Union{Vector{spg_encodetype}, Vector{spg_symtype}}
    encodeds = [nspg_symmetry_operations[i] for i = nspg_index[nspgidx]:nspg_index[nspgidx]+nspg_data[nspgidx].order-1]
    if encode
        return encodeds
    else
        syms = [spg_decode_symmetry(encoded) for encoded in encodeds]
        return syms
    end
end

"""
get_full_symmetry_operations_of_spin_point_group(spgidx::Int64; encode::Bool) -> Union{Vector{spg_encodetype}, Vector{spg_symtype}}
given a spin point group number, find all its symmetry operations
"""
function get_full_symmetry_operations_of_spin_point_group(spgidx::Int64; encode::Bool)::Union{Vector{spg_encodetype}, Vector{spg_symtype}}
    encodeds = [spg_symmetry_operations[i] for i = spg_index[spgidx]:spg_index[spgidx]+spg_data[spgidx].order-1]
    if encode
        return encodeds
    else
        syms = [spg_decode_symmetry(encoded) for encoded in encodeds]
        return syms
    end
end

function get_reduced_group(encodeds::Vector{spg_encodetype}; SO2::Int64, encode::Bool)
    syms = [spg_decode_symmetry(encoded) for encoded in encodeds]
    return get_reduced_group(syms, encode=encode, SO2=SO2)
end
function get_reduced_group(syms::Vector{spg_symtype}; SO2::Int64, encode::Bool)
    nsyms = get_reduced_group(syms, SO2)
    if encode
        nencodeds = [spg_encode_symmetry(nsym) for nsym in nsyms]
        return nencodeds
    else
        return nsyms
    end
end
function get_reduced_group(syms::Vector{spg_symtype}, SO2::Int64)::Vector{spg_magnonbasis_symtype}
    G = spg_magnonbasis_symtype[
        (r=[1 0 0; 0 1 0; 0 0 1], s=[ 1 0; 0  1], T=1),
        (r=[1 0 0; 0 1 0; 0 0 1], s=[ 1 0; 0 -1], T=1),
        (r=[1 0 0; 0 1 0; 0 0 1], s=[-1 0; 0 -1], T=1),
        (r=[1 0 0; 0 1 0; 0 0 1], s=[-1 0; 0  1], T=1)
    ]
    reducedgroup = spg_magnonbasis_symtype[]
    nsyms = spg_magnonbasis_symtype[change_symmetry_from_spin_basis_to_magnon_basis(sym, SO2=SO2) for sym in syms]
    for g in G
        for sym in nsyms
            result = composite_symmetries(g, sym, encode=false)
            push!(reducedgroup, result)
        end
    end
    if !isa_group(reducedgroup) error("not a group") end
    return reducedgroup
end

"""
isa_set(list::Vector{<:Union{spg_symtype, spg_encodetype, spg_magnonbasis_symtype, spg_magnonbasis_encodetype}}) -> Bool
check if a list forms a set,i.e.,each element only appears once
"""
function isa_set(list::Vector{<:Union{spg_symtype, spg_encodetype, spg_magnonbasis_symtype, spg_magnonbasis_encodetype}})::Bool
    repeated=[]
    for i in eachindex(list)
        for j in eachindex(list)
            if i == j continue end
            if is_symmetry_same(list[i], list[j]) return false end
        end
    end
    return true
end

"""
isa_group(encodeds::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}}) -> Bool
isa_group(syms::Vector{spg_symtype}) -> Bool
isa_group(syms::Vector{spg_magnonbasis_symtype}) -> Bool
check if a set forms a group
"""
function isa_group(encodeds::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}})::Bool
    syms = [spg_decode_symmetry(encoded) for encoded in encodeds]
    return isa_group(syms)
end
function isa_group(syms::Vector{spg_symtype})::Bool
    if !isa_set(syms) return false end
    return sptsymmetry.isa_group([(r=sym.r, s=Float64.(sym.s)) for sym in syms])
end
function isa_group(syms::Vector{spg_magnonbasis_symtype})::Bool
    if !isa_set(syms) return false end
    return sptsymmetry.isa_group([(r=sym.r, s=Float64.(sym.s), T=sym.T) for sym in syms])
end

"""
generate_full_spin_point_group_of_nontrivial_spin_point_group(nspgidx, arrgmt)
"""
function generate_full_spin_point_group_of_nontrivial_spin_point_group(nspgidx, arrgmt)
    nspgsyms = get_full_symmetry_operations_of_nontrivial_spin_point_group(nspgidx, encode=false)
    nspgnot = get_notation_of_nontrivial_spin_point_group(nspgidx, subscript=true)
    holo = split(string(nspg_data[nspgidx].Bholohedry), ".")[end]
    if !(holo in ["TRIGO", "HEXA"])
        rms = [
            [-1 0 0;0 1 0;0 0 1],
            [1 0 0;0 -1 0;0 0 1],
            [1 0 0;0 1 0;0 0 -1]
        ]
    else
        rms = [
            [-1 1 0;0 1 0;0 0 1],
            [1 0 0;1 -1 0;0 0 1],
            [1 0 0;0 1 0;0 0 -1]
        ]
    end
    r6s = [
        [1 0 0; 0 1 -1; 0 1 0],
        [-1 0 1; 0 1 0; -1 0 0],
        [1 -1 0; 1 0 0; 0 0 1]
    ]
    spginfos = []
    for result in arrgmt
        grouptype = result[1]
        spintype = result[2]
        constraint = result[3]
        SO2dir = nothing
        if grouptype != 1
            if !haskey(constraint, "relation")
                if spintype == [3]
                    sogm = "mz"
                else
                    sogm = "mx"
                    sogr = "z"
                end
            elseif constraint["relation"] == "∥"
                if [1,0,0] == constraint["axis"]
                    if spintype == [3]
                        sogm = "mx"
                    else
                        SO2dir = "∥"
                        sogm = "mz"
                        sogr = "x"
                    end
                elseif [0,0,1] == constraint["axis"]
                    if spintype == [3]
                        sogm = "mz"
                    else
                        SO2dir = "∥"
                        sogm = "mx"
                        sogr = "z"
                    end
                elseif [0,1,0] == constraint["axis"]
                    if spintype == [3]
                        sogm = "my"
                    else
                        SO2dir = "∥"
                        sogm = "mz"
                        sogr = "y"
                    end
                end
            elseif constraint["relation"] == "⟂"
                if [0,0,1] == constraint["axis"]
                    if spintype == [3]
                        sogm = "mx"
                    else
                        SO2dir = "⟂"
                        sogm = "mz"
                        sogr = "x"
                    end
                elseif [0,1,0] == constraint["axis"]
                    if spintype == [3]
                        sogm = "mz"
                    else
                        SO2dir = "⟂"
                        sogm = "mz"
                        sogr = "x"
                    end
                else
                    error(123)
                end
            else
                error(123)
            end
        end
        if grouptype == 1
            mats = deepcopy(nspgsyms)
            spginfo = (notation = nspgnot, operations = mats, order = length(mats), grouptype=1, nspg=nspgidx, SO2=0, SO2dir=SO2dir, C2T=0, spinconfig=spintype)
        elseif grouptype == 2
            if sogm == "mz"
                rm = (r=[1 0 0;0 1 0;0 0 1], s=rms[3])
                midx = 3
            elseif sogm == "mx"
                rm = (r=[1 0 0;0 1 0;0 0 1], s=rms[1])
                midx = 1
            elseif sogm == "my"
                rm = (r=[1 0 0;0 1 0;0 0 1], s=rms[2])
                midx = 2
            end
            mats = [nspgsyms; [composite_symmetries(rm, sym, encode=false) for sym in nspgsyms]]
            not = "$(nspgnot)($(sogm))1"
            spginfo = (notation = not, operations = mats, order = length(mats), grouptype=2, nspg=nspgidx, SO2=0, SO2dir=SO2dir, C2T=midx, spinconfig=spintype)
            #=
            isotable = get_isometry_type_table(spginfo.operations)
            for otherspginfo in spginfos
                if otherspginfo.grouptype != 2 continue end
                otherisotable = get_isometry_type_table(otherspginfo.operations)
                if all(isotable.==otherisotable) @goto nextspginfo end
            end
            =#
        elseif grouptype == 3
            mats = deepcopy(nspgsyms)
            not = "$(nspgnot)(∞$(sogr))1"
            if sogr == "x"
                so2idx = 1
            elseif sogr == "y"
                so2idx = 2
            elseif sogr == "z"
                so2idx = 3
            end
            spginfo = (notation = not, operations = mats, order = length(mats), grouptype=3, nspg=nspgidx, SO2=so2idx, SO2dir=SO2dir, C2T=0, spinconfig=spintype)
            #=
            isotable = get_isometry_type_table(get_reduced_operations_of_spin_point_group_with_SO2(spginfo.operations, spginfo.SO2))
            for otherspginfo in spginfos
                if otherspginfo.grouptype != 3 continue end
                otherisotable = get_isometry_type_table(get_reduced_operations_of_spin_point_group_with_SO2(otherspginfo.operations, otherspginfo.SO2))
                if all(isotable.==otherisotable) @goto nextspginfo end
            end
            =#
        elseif grouptype == 4
            if sogm == "mz"
                rm = (r=[1 0 0;0 1 0;0 0 1], s=rms[3])
                midx = 3
            elseif sogm == "mx"
                rm = (r=[1 0 0;0 1 0;0 0 1], s=rms[1])
                midx = 1
            elseif sogm == "my"
                rm = (r=[1 0 0;0 1 0;0 0 1], s=rms[2])
                midx = 2
            end
            mats = [nspgsyms; [composite_symmetries(rm, sym, encode=false) for sym in nspgsyms]]
            not = "$(nspgnot)(∞$(sogr)$(sogm))1"
            if sogr == "x"
                so2idx = 1
            elseif sogr == "y"
                so2idx = 2
            elseif sogr == "z"
                so2idx = 3
            end
            spginfo = (notation = not, operations = mats, order = length(mats), grouptype=4, nspg=nspgidx, SO2=so2idx, SO2dir=SO2dir, C2T=midx, spinconfig=spintype)
        end
        push!(spginfos, spginfo)
        @label nextspginfo
    end
    return spginfos
end

end #(module GroupInformation)

#################################################GroupRepresentation##################################################
######################################################################################################################
module GroupRepresentation
using LinearAlgebra
using ...MathTool
using ...GroupRepresentation
using ...PointGroup.GroupRepresentation:get_character_table as get_character_table_of_point_group
using ...PointGroup.GroupRepresentation:get_irreps as get_irreps_of_point_group
using ..SpinPointGroupDataBase
using ..SpinPointSymmetry
using ..GroupInformation

export get_character_table_of_group
export get_irreps_of_group
export get_irreps_of_direct_product
export get_induced_irreps
export get_character_table_of_finite_full_spin_point_group
export get_irreps_of_finite_full_spin_point_group
export get_character_table_of_nontrivial_spin_point_group
export get_irreps_of_nontrivial_spin_point_group


function get_character_table_of_group(syms::Vector{<:Union{spg_symtype, spg_magnonbasis_symtype}})
    encodeds = [spg_decode_symmetry(sym) for sym in syms]
    return get_character_table_of_group(encodeds)
end
function get_character_table_of_group(encodeds::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}})
    irreps, case = get_irreps_of_group(encodeds, coirrep=false)
    character_table, classes = get_character_table_of_irreps(irreps, class=true)
    return (table=character_table, class=classes, case=case)

end

"""
get_irreps_of_group(syms::Vector{spg_magnonbasis_symtype}; coirrep::Bool)
get_irreps_of_group(encodeds::Vector{spg_magnonbasis_encodetype}; coirrep::Bool)
formula
    nspg = P ⊕ qP

    fspg
    = nspg ⊗ {[E||E], [C⟂T||E]}
    = (P ⊕ qP) ⊗ {[E||E], [C⟂T||E]}
    fspg
    = U ⊕ aU
    # U: unitary subgroup of fspg, a: antiunitary sym of fspg
    = (P' ⊕ qP') ⊕ a(P' ⊕ qP')

    colAFM
    = fspg  ⊗ {[E||E], [Σz||E], [-E||E], [-Σz||E]}
    = ((P' ⊕ qP') ⊕ a(P' ⊕ qP')) ⊗ {[E||E], [Σz||E], [-E||E], [-Σz||E]}
    colAFM
    = Ugroup ⊕ aUgroup
    = (Pgroup ⊕ Qgroup) ⊕ a(Pgroup ⊕ Qgroup)

    Pgroup = P' ⊗ {[E||E], [Σz||E], [-E||E], [-Σz||E]}
    Qgroup = qP' ⊗ {[E||E], [Σz||E], [-E||E], [-Σz||E]}
    Ugroup = (P' ⊕ qP') ⊗ {[E||E], [Σz||E], [-E||E], [-Σz||E]}

procedure
    1.get irrep of Pgroup
    2.get induced irrep and character table of Ugroup = Pgroup ⊕ Qgroup
    3.get cocharacter table of colAFM = Ugroup ⊕ aUgroup
"""
function get_irreps_of_group(syms::Vector{spg_symtype}; coirrep::Bool)::NamedTuple{(:irreps, :case), Tuple{Matrix{Matrix{ComplexF64}}, Vector{Int64}}}
    encodeds =[spg_encode_symmetry(sym) for sym in syms]
    return get_irreps_of_group(encodeds, coirrep=coirrep)
end
function get_irreps_of_group(encodeds::Vector{spg_encodetype}; coirrep::Bool)::NamedTuple{(:irreps, :case), Tuple{Matrix{Matrix{ComplexF64}}, Vector{Int64}}}
    Ugroup = get_maximal_unitary_subgroup(encodeds, encode = true)
    aUgroup = get_maximal_antiunitary_set(encodeds, encode = true)
    if aUgroup == []
        #=get pg isomorphic to nspg=#
        ptgroup = Vector{Int64}([encoded[1] for encoded in encodeds])
        #=get irreps of pg/nspg=#
        irreps = get_irreps_of_point_group(ptgroup)
        if !check_irreps(irreps) error("error occurs at get_irreps_of_group: irreps not correct") end
        return (irreps=irreps, case=Int64[])
    else
        irreps, _ = get_irreps_of_group(Ugroup, coirrep=coirrep)
        return get_coirreps(encodeds, Ugroup, aUgroup, irreps, coirrep=coirrep)
    end
end
function get_irreps_of_group(syms::Vector{spg_magnonbasis_symtype}; coirrep::Bool)
    encodeds = [spg_encode_symmetry(sym) for sym in syms]
    return get_irreps_of_group(encodeds, coirrep=coirrep)
end
function get_irreps_of_group(encodeds::Vector{spg_magnonbasis_encodetype}; coirrep::Bool)
    #= get Ugroup, aUgroup, Pgroup and Qgroup=#
    Ugroup = get_maximal_unitary_subgroup(encodeds, encode=true)
    aUgroup = get_maximal_antiunitary_set(encodeds, encode=true)
    Pgroup = spg_magnonbasis_encodetype[]
    Qgroup = spg_magnonbasis_encodetype[]
    for encoded in Ugroup
        sym = spg_decode_symmetry(encoded)
        if all(Diagonal(sym.s) .== 0)
            push!(Qgroup, encoded)
        else
            push!(Pgroup, encoded)
        end
    end
    #=get irrep of Pgroup=#
    irreps_list = get_irreps_of_direct_product(Pgroup)
    #=get induced irrep of Ugroup=#
    if length(Qgroup) != 0 irreps_list = get_induced_irreps(Ugroup, Pgroup, irreps_list) end
    #==get coirrep=#
    if length(aUgroup) != 0
        irreps_list, case = get_coirreps(encodeds, Ugroup, aUgroup, irreps_list, coirrep=coirrep)
    else
        case = Int64[]
    end

    #=remove unphysical irrep=#
    idx = findall(x->x==(16484, 12, 1), encodeds)[1] #\bar{E}
    idxs = [i for i in axes(irreps_list, 1) if all(real(diag(irreps_list[i, idx])) .< 0.0)]
    irreps = Matrix{ComplexF64}[irreps_list[i, j] for i in idxs, j in axes(irreps_list, 2)]
    if length(aUgroup) != 0 && coirrep case = [case[i] for i in idxs] end

    return (irreps=irreps, case=case)

end

"""
formula
    group = nspg ⊗ {[E||E], [Σz||E], [-E||E], [-Σz||E]}
procedure
    1.get pg isomorphic to nspg
    2.get irreps of pg/nspg
    3.get irreps of direct product of nspg and {[E||E], [Σz||E], [-E||E], [-Σz||E]}
directproduct_irreps
        encoded[1] encoded[2] encoded[3] ...
irrep1
irrep2
...
"""
function get_irreps_of_direct_product(encodeds::Vector{spg_magnonbasis_encodetype})::Matrix{Matrix{Complex{Float64}}}
    pointgroup = Vector{Int64}()
    #=get pg isomorphic to nspg=#
    for encoded in encodeds
        if !(encoded[1] in pointgroup) push!(pointgroup, encoded[1]) end
    end
    #=get irreps of pg/nspg=#
    irreps = get_irreps_of_point_group(pointgroup)
    #=get irreps of direct product of nspg and {[E||E], [Σz||E], [-E||E], [-Σz||E]}=#
    directproduct_irreps = Matrix{Any}([nothing for i=1:4*size(irreps, 1), j=1:4*size(irreps, 2)])
    for j in eachindex(encodeds)
        encoded = encodeds[j]
        for i in eachindex(pointgroup)
            if encoded[1] == pointgroup[i]
                if encoded[2] == 68
                    directproduct_irreps[:, j] = repeat(irreps[:, i], 4)
                elseif encoded[2] == 66
                    directproduct_irreps[:, j] = [irreps[:, i]; irreps[:, i].*(-1); irreps[:, i]; irreps[:, i].*(-1)]
                elseif encoded[2] == 12
                    directproduct_irreps[:, j] = [irreps[:, i]; irreps[:, i]; irreps[:, i].*(-1); irreps[:, i].*(-1)]
                elseif encoded[2] == 14
                    directproduct_irreps[:, j] = [irreps[:, i]; irreps[:, i].*(-1); irreps[:, i].*(-1); irreps[:, i]]
                else
                    error("Error occurs at get_irreps")
                end
                break
            end
        end
    end
    for i in axes(directproduct_irreps, 1), j in axes(directproduct_irreps, 2)
        if size(directproduct_irreps[i, j]) == () directproduct_irreps[i, j] = [directproduct_irreps[i, j];;] end
    end
    directproduct_irreps = Matrix{Matrix{Complex{Float64}}}(directproduct_irreps)
    if !check_irreps(directproduct_irreps) error("irrep not correct") end
    return directproduct_irreps
end

"""
get induced irrep of super from sub(only solve invariant subgroup of index 2)
formula
    super = sub + q*sub
    D(h)_i = ith irrep of sub
    (D(h)_i)q = D(q^-1*h*q)_i = conjugate irrep ~ D(h)_j = jth irrep of sub

    if i == j or (D(h)_i)q ~ D(h)_i
        D'(h) = D(h)_i
        D'(q) = U
        U^-1 * D(h)_i * U = (D(h)_i)q
    if i != j or (D(h)_i)q ~ D(h)_j
        D'(h) = [
            D(h)_i      0
              0     (D(h)_i)q
        ]
        D'(q) = [
            0 D(q^2)_i
            I    0
        ]
procedure
    1. get conjugate irrepresntation for each irrepresntation
    2. induce irrep depending on if irrep and conj_irrep are similar
"""
function get_induced_irreps(super::Vector{spg_magnonbasis_encodetype}, sub::Vector{spg_magnonbasis_encodetype}, irreps_list::Matrix{Matrix{Complex{Float64}}}; tor=1e-3)::Matrix{Matrix{Complex{Float64}}}
    irreps_list = convert(Matrix{Matrix{Complex{Float64}}}, irreps_list)
    if length(super) != 2* length(sub) error("only solve invariant subgroup of index 2") end

    q = [encoded for encoded in super if !(encoded in sub)][1]#coset representative
    conj_sub = [conjugate_symmetry(encoded, q, encode=true) for encoded in sub]
    conj_idx = [findall(x->x==encoded, sub)[1] for encoded in conj_sub]
    gorder = size(irreps_list, 2)
    used_irrep_idx = []
    ind_irreps_list = []
    for i in axes(irreps_list, 1)
        if i in used_irrep_idx continue end
        #=get conjugate irrepresentation=#
        conj_irreps = [irreps_list[i, idx] for idx in conj_idx]
        #=induce irrep depending on if irrep and conj_irrep are similar=#
        j = [k for k in axes(irreps_list, 1) if all(abs.(tr.(conj_irreps) - tr.(irreps_list[k, :])) .< tor)][1]
        dim = size(irreps_list[i, 1], 1)
        if i == j
            Tmatrix = get_similar_transformation_matrix_between_two_reps(irreps_list[i, :], conj_irreps)
            q_irreps = [Tmatrix, -Tmatrix]
            ind_irreps = Matrix{Matrix{Complex{Float64}}}([zeros(dim, dim) for i in 1:2, j in 1:length(super)])
            for pmidx = 1:2
                q_irrep = q_irreps[pmidx]
                for k in eachindex(sub)
                    u = sub[k]
                    u_irrep = Matrix{Complex{Float64}}(irreps_list[i, k])
                    Idx = findall(x->x==u, super)[1]
                    ind_irreps[pmidx, Idx] = u_irrep
                    qu = composite_symmetries(q, u, encode=true)
                    qu_irrep = composite_matrice(q_irrep, u_irrep)
                    Idx = findall(x->x==qu, super)[1]
                    ind_irreps[pmidx, Idx] = qu_irrep
                end
            end
        else
            q2_idx = findall(x->x==composite_symmetries(q, order=2, encode=true), sub)[1]
            q_irrep = Matrix{Complex{Float64}}([[zeros(dim, dim); Matrix(I, dim, dim)];; [irreps_list[i, q2_idx]; zeros(dim, dim)]])
            ind_irreps = Matrix{Matrix{Complex{Float64}}}([zeros(2*dim, 2*dim) for i in 1:1, j in 1:length(super)])
            for k in eachindex(sub)
                u = sub[k]
                u_irrep = Matrix{Complex{Float64}}([[irreps_list[i, k]; zeros(dim, dim)];; [zeros(dim, dim); conj_irreps[k]]])
                Idx = findall(x->x==u, super)[1]
                ind_irreps[1, Idx] = u_irrep
                qu = composite_symmetries(q, u, encode=true)
                qu_irrep = composite_matrice(q_irrep, u_irrep)
                Idx = findall(x->x==qu, super)[1]
                ind_irreps[1, Idx] = qu_irrep
            end
            append!(used_irrep_idx,[i, j])
            @label next
        end
        if ind_irreps_list == []
            ind_irreps_list = ind_irreps
        else
            ind_irreps_list = [ind_irreps_list; ind_irreps]
        end
    end
    return ind_irreps_list
end

"""
#=non-unitary group, get coirrepresentation=#
"""
function get_coirreps(encodeds::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}}, Ugroup::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}}, aUgroup::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}}, irreps::Matrix{Matrix{ComplexF64}}; coirrep::Bool)#::NamedTuple{(:irreps, :case), Tuple{Matrix{Matrix{ComplexF64}}, Vector{Int64}}}#Tuple{Matrix{Matrix{ComplexF64}}, Vector{Int64}}
    #=get Wigner criterion=#
    character_sum = zeros(ComplexF64, size(irreps, 1))
    for sym in aUgroup
        sym_squared = composite_symmetries(sym, order=2, encode = true)
        idx = findall(x->x==sym_squared, Ugroup)[1]
        for i = 1:size(irreps, 1)
            character_sum[i] += tr(irreps[i, idx])/length(Ugroup)
        end
    end
    character_sum = round.(character_sum, digits=3)
    if all(isinteger.(character_sum))
        character_sum = Int64.(character_sum)
    else
        error("Error occurs at get_irreps_of_nontrivial_spin_point_group: character sum is not a real integer")
    end

    if !coirrep
        case = []
        a0 = aUgroup[1]
        a0_inv = get_inverse_of_symmetry(a0, encode=true)
        used_idx = []
        for i in 1:size(irreps, 1)
            if i in used_idx continue end
            conj_irrep = [conj(irreps[i, findall(x->x==composite_symmetries(a0_inv, u, a0, encode=true), Ugroup)[1]]) for u in Ugroup]
            if character_sum[i] == 1
                push!(case, (1, i))
            elseif character_sum[i] == 0
                j = get_index_of_irrep(conj_irrep, irreps)
                push!(used_idx, j)
                push!(case, (0, i, j))
            elseif character_sum[i] == -1
                push!(case, (-1, i))
            end
        end
        return (irreps=irreps, case=case)
    else
        #=get coirreps=#
        case = Int64[]
        coirreps = Matrix{ComplexF64}[]
        a0 = aUgroup[1]
        a0_inv = get_inverse_of_symmetry(a0, encode=true)
        used_idx = []
        for i in 1:size(irreps, 1)
            if i in used_idx continue end
            conj_irrep = [conj(irreps[i, findall(x->x==composite_symmetries(a0_inv, u, a0, encode=true), Ugroup)[1]]) for u in Ugroup]
            if size(irreps[i, 1]) == ()
                zeromat = ComplexF64(0)
            else
                zeromat = zeros(ComplexF64, size(irreps[i, 1]))
            end
            coirrepi = Matrix{ComplexF64}[]
            if character_sum[i] == 1
                a0mat = get_similar_transformation_matrix_between_two_reps(irreps[i,:], conj_irrep)
                for encoded in encodeds
                    if is_symmetry_unitary(encoded)
                        idx = findall(x->x==encoded, Ugroup)[1]
                        push!(coirrepi, irreps[i, idx])
                    else
                        a = encoded
                        idx = findall(x->x==composite_symmetries(a, a0_inv, encode=true), Ugroup)[1]
                        amat = irreps[i, idx] * a0mat
                        push!(coirrepi, amat)
                    end
                end
            elseif character_sum[i] == 0
                push!(used_idx, get_index_of_irrep(conj_irrep, irreps))
                for encoded in encodeds
                    if is_symmetry_unitary(encoded)
                        u = encoded
                        idx1 = findall(x->x==encoded, Ugroup)[1]
                        idx2 = findall(x->x==composite_symmetries(a0_inv, u, a0, encode=true), Ugroup)[1]
                        umat = [[irreps[i, idx1] ;; zeromat] ; [zeromat ;; conj(irreps[i, idx2])]]
                        push!(coirrepi, umat)
                    else
                        a = encoded
                        idx1 = findall(x->x==composite_symmetries(a, a0, encode=true), Ugroup)[1]
                        idx2 = findall(x->x==composite_symmetries(a0_inv, a, encode=true), Ugroup)[1]
                        amat = [[zeromat ;; irreps[i, idx1]] ; [conj(irreps[i, idx2]) ;; zeromat]]
                        push!(coirrepi, amat)
                    end
                end
            elseif character_sum[i] == -1
                a0mat = get_similar_transformation_matrix_between_two_reps(irreps[i,:], conj_irrep)
                for encoded in encodeds
                    if is_symmetry_unitary(encoded)
                        idx = findall(x->x==encoded, Ugroup)[1]
                        umat = [[irreps[i, idx] ;; zeromat] ; [zeromat ;; irreps[i, idx]]]
                        push!(coirrepi, umat)
                    else
                        a = encoded
                        idx = findall(x->x==composite_symmetries(a, a0_inv, encode=true), Ugroup)[1]
                        mat = irreps[i, idx]*a0mat
                        amat = [[zeromat ;; mat] ; [mat ;; zeromat]]
                        push!(coirrepi, amat)
                    end
                end
            end
            push!(case, character_sum[i])
            if coirreps == []
                coirreps = reshape(coirrepi, 1, length(coirrepi))
            else
                coirreps = [coirreps ; reshape(coirrepi, 1, length(coirrepi))]
            end
        end
        return (irreps=coirreps, case=case)
    end
end

end #(module GroupRepresentation)

########################################################MagneticMoment########################################################
######################################################################################################################
module MagneticMoment
using LinearAlgebra
using ..SpinPointGroupDataBase
using ..GroupInformation
import ...PointSymmetry as ptsym
using ...PointGroup.GroupInformation:get_class_of_point_group

export get_supporting_magnetic_system
export get_supporting_spin_arrangement_of_nosoc

function get_supporting_magnetic_system(nspgidx; tor=1e-3)
    syms = get_full_symmetry_operations_of_nontrivial_spin_point_group(nspgidx, encode=false)
    main_axis = (6, nothing)
    order2_axis = Dict(4=>Vector{Tuple{Vector{Int64}, Vector{Int64}}}(), 7=>Vector{Tuple{Vector{Int64}, Vector{Int64}}}())
    ssyms = Vector{Matrix{Int64}}()
    for sym in syms
        s = sym.s
        if !(s in ssyms) push!(ssyms, s) end
        iso_type = ptsym.get_isometry_type_of_point_symmetry(s)
        iso_axis = ptsym.get_isometry_axis_of_point_symmetry(s)
        if iso_type == 6 continue end
        if iso_type in [4, 7]
            r_iso_type = ptsym.get_isometry_type_of_point_symmetry(sym.r)
            for i in eachindex(order2_axis[iso_type])
                (axis, r_iso_types) = order2_axis[iso_type][i]
                if (iso_axis == axis || -iso_axis == axis)
                    if !(r_iso_type in r_iso_types) push!(order2_axis[iso_type][i][2], r_iso_type) end
                    @goto next
                end
            end
            push!(order2_axis[iso_type], (iso_axis, [r_iso_type]))
        end
        @label next
        iso_type = ptsym.get_isometry_type_of_point_symmetry(Int64.(det(s)*s))
        if iso_type-6 > main_axis[1]-6 main_axis = (iso_type, iso_axis) end
    end
    pgidx = get_class_of_point_group(ssyms)
    if pgidx >= 28
        #multiple high order main axis
        arrgmt = [
            (1, [4], Dict()) #Gnspg
        ]
    elseif pgidx in [9, 16, 21]
        #one high order main axis, proper group
        arrgmt = [
            (1, [3, 4], Dict()), #Gnspg
            (2, [3], Dict("axis"=>main_axis[2], "relation"=>"∥")), #Gnspg x Z^{k}_2
            (3, [1, 2], Dict("axis"=>main_axis[2], "relation"=>"∥")), #Gnspg x ∞z
            (4, [1, 2], Dict("axis"=>main_axis[2], "relation"=>"∥")) #Gnspg x {∞zm⟂}1
        ]
    elseif pgidx in [10, 11, 17, 22, 23]
        #one high order main axis, improper group
        arrgmt = [
            (1, [3, 4], Dict()), #Gnspg
            (2, [3], Dict("axis"=>main_axis[2], "relation"=>"∥")), #Gnspg x Z^{k}_2
            (3, [2], Dict("axis"=>main_axis[2], "relation"=>"∥")), #Gnspg x ∞z
            (4, [2], Dict("axis"=>main_axis[2], "relation"=>"∥")) #Gnspg x {∞zm⟂}1
        ]
    elseif pgidx in [12, 13, 14, 15, 18, 19, 20, 24, 25, 26, 27]
        #one high order main axis, multiple twofold axis
        arrgmt = [
            (1, [3, 4], Dict()), #Gnspg
            (2, [3], Dict("axis"=>main_axis[2], "relation"=>"∥")), #Gnspg x Z^{k}_2
            (3, [2], Dict("axis"=>main_axis[2], "relation"=>"∥")), #Gnspg x ∞z
            (4, [2], Dict("axis"=>main_axis[2], "relation"=>"∥")) #Gnspg x {∞zm⟂}1
        ]
    elseif pgidx in [6, 7, 8]
        #multiple twofold axis
        arrgmt = []
        #Gnspg
        push!(arrgmt, (1, [3, 4], Dict()))
        #Gnspg x Z^{k}_2
        for i in keys(order2_axis)
            j = 1
            while j <= length(order2_axis[i])
                for k = 1:j-1
                    if sort!(order2_axis[i][j][2]) == sort!(order2_axis[i][k][2])
                        deleteat!(order2_axis[i], j)
                        @goto next_j
                    end
                end
                j += 1
                @label next_j
            end
        end
        arrgmt = [arrgmt;
            [
            (2, [3], Dict("axis"=>axis, "relation"=>"∥"))
            for axes in values(order2_axis)
            for (axis, r_iso_type) in axes
            ]
        ]
        #Gnspg x ∞z
        arrgmt = [arrgmt;
            [
            (3, [2], Dict("axis"=>axis, "relation"=>"∥"))
            for axes in values(order2_axis)
            for (axis, r_iso_type) in axes
            ]
        ]
        #Gnspg x {∞zm⟂}1
        arrgmt = [arrgmt;
            [
            (4, [2], Dict("axis"=>axis, "relation"=>"∥"))
            for axes in values(order2_axis)
            for (axis, r_iso_type) in axes
            ]
        ]
    elseif pgidx == 4
        #m
        arrgmt = [
            (1, [3, 4], Dict()), #Gnspg
            (2, [3], Dict("axis"=>main_axis[2], "relation"=>"∥")), #Gnspg x Z^{k}_2
            (2, [3], Dict("axis"=>main_axis[2], "relation"=>"⟂")), #Gnspg x Z^{k}_2
            (3, [2], Dict("axis"=>main_axis[2], "relation"=>"∥")), #Gnspg x ∞z
            (3, [1, 2], Dict("axis"=>main_axis[2], "relation"=>"⟂")), #Gnspg x ∞z
            (4, [2], Dict("axis"=>main_axis[2], "relation"=>"∥")), #Gnspg x {∞zm⟂}1
            (4, [1, 2], Dict("axis"=>main_axis[2], "relation"=>"⟂")) #Gnspg x {∞zm⟂}1
        ]
    elseif pgidx == 5
        #2/m
        arrgmt = [
            (1, [3, 4], Dict()), #Gnspg
            (2, [3], Dict("axis"=>main_axis[2], "relation"=>"∥")), #Gnspg x Z^{k}_2
            (2, [3], Dict("axis"=>main_axis[2], "relation"=>"⟂")), #Gnspg x Z^{k}_2
            (3, [2], Dict("axis"=>main_axis[2], "relation"=>"∥")), #Gnspg x ∞z
            (3, [2], Dict("axis"=>main_axis[2], "relation"=>"⟂")), #Gnspg x ∞z
            (4, [2], Dict("axis"=>main_axis[2], "relation"=>"∥")), #Gnspg x {∞zm⟂}1
            (4, [2], Dict("axis"=>main_axis[2], "relation"=>"⟂")) #Gnspg x {∞zm⟂}1
        ]
     elseif pgidx == 2
        #-1
        arrgmt = [
            (1, [3, 4], Dict()), #Gnspg
            (2, [3], Dict("axis"=>[0, 0, 0])), #Gnspg x Z^{k}_2
            (3, [2], Dict("axis"=>[0, 0, 0])), #Gnspg x ∞z
            (4, [2], Dict("axis"=>[0, 0, 0])) #Gnspg x {∞zm⟂}1
        ]
    elseif pgidx == 3
        #2
        arrgmt = [
            (1, [3, 4], Dict()), #Gnspg
            (2, [3], Dict("axis"=>main_axis[2], "relation"=>"∥")), #Gnspg x Z^{k}_2
            (2, [3], Dict("axis"=>main_axis[2], "relation"=>"⟂")), #Gnspg x Z^{k}_2
            (3, [1, 2], Dict("axis"=>main_axis[2], "relation"=>"∥")), #Gnspg x ∞z
            (3, [2], Dict("axis"=>main_axis[2], "relation"=>"⟂")), #Gnspg x ∞z
            (4, [1, 2], Dict("axis"=>main_axis[2], "relation"=>"∥")), #Gnspg x {∞zm⟂}1
            (4, [2], Dict("axis"=>main_axis[2], "relation"=>"⟂")) #Gnspg x {∞zm⟂}1
        ]
    elseif pgidx == 1
        #1
        arrgmt = [
            (1, [3, 4], Dict()), #Gnspg
            (2, [3], Dict("axis"=>[0, 0, 0])), #Gnspg x Z^{k}_2
            (3, [1, 2], Dict("axis"=>[0, 0, 0])), #Gnspg x ∞z
            (4, [1, 2], Dict("axis"=>[0, 0, 0])) #Gnspg x {∞zm⟂}1
        ]
    end
    return arrgmt
end

"""
get_supporting_spin_arrangement_of_nosoc(idx::Int64) -> Vector{Tuple{Int64, Dict{String, Any}}}
get_supporting_spin_arrangement_of_nosoc(encodeds::Vector{spg_encodetype}) -> Vector{Tuple{Int64, Dict{String, Any}}}
get_supporting_spin_arrangement_of_nosoc(syms::Vector{spg_symtype}; tor=1e-3) -> Vector{Tuple{Int64, Dict{String, Any}}}

Return: arrgmt
    arrgmt[i] = (type, constraint)
    type=1/2 : collinearFM/collinearAFM
        constraint=Dict("axis"=>[0, 0, 0]) : no constraint on m
        constraint=Dict("axis"=>vec, "relation"=>"∥") : m ∥ vec
        *m:direction vector of spin line, i.e., spin direction ∥ m
    type=3 : coplanar
        constraint=Dict("axis"=>vec, "relation"=>"∥") : n ∥ vec
        constraint=Dict("axis"=>vec, "relation"=>"⟂") : n ⟂ vec
        *n:normal vector of spin plane, i.e., spin direction ⟂ n
    type=4 : non-coplanar
        constraint=Dict() : constraint depends on specific group
"""
function get_supporting_spin_arrangement_of_nosoc(idx::Int64)::Vector{Tuple{Int64, Dict{String, Any}}}
    syms = get_full_symmetry_operations_of_nontrivial_spin_point_group(idx, encode=false)
    return get_supporting_spin_arrangement_of_nosoc(syms)
end
function get_supporting_spin_arrangement_of_nosoc(encodeds::Vector{spg_encodetype})::Vector{Tuple{Int64, Dict{String, Any}}}
    syms = [spg_decode_symmetry(encoded) for encoded in encodeds]
    return get_supporting_spin_arrangement_of_nosoc(syms)
end
function get_supporting_spin_arrangement_of_nosoc(syms::Vector{spg_symtype}; tor=1e-3)::Vector{Tuple{Int64, Dict{String, Any}}}
    main_axis = (6, nothing)
    order2_axis = Dict(4=>Vector{Tuple{Vector{Int64}, Vector{Int64}}}(), 7=>Vector{Tuple{Vector{Int64}, Vector{Int64}}}())
    ssyms = Vector{Matrix{Int64}}()
    for sym in syms
        s = sym.s
        if !(s in ssyms) push!(ssyms, s) end
        iso_type = get_isometry_type_of_point_symmetry(s)
        iso_axis = get_isometry_axis_of_point_symmetry(s)
        if iso_type == 6 continue end
        if iso_type in [4, 7]
            r_iso_type = get_isometry_type_of_point_symmetry(sym.r)
            for i in eachindex(order2_axis[iso_type])
                (axis, r_iso_types) = order2_axis[iso_type][i]
                if (iso_axis == axis || -iso_axis == axis)
                    if !(r_iso_type in r_iso_types) push!(order2_axis[iso_type][i][2], r_iso_type) end
                    @goto next
                end
            end
            push!(order2_axis[iso_type], (iso_axis, [r_iso_type]))
        end
        @label next
        iso_type = get_isometry_type_of_point_symmetry(Int64.(det(s)*s))
        if iso_type-6 > main_axis[1]-6 main_axis = (iso_type, iso_axis) end
    end
    pgidx = get_class_of_point_group(ssyms)
    if pgidx >= 28
        #multiple high order main axis
        arrgmt = [
            (4, Dict())
        ]
    elseif pgidx in [12, 13, 14, 15, 18, 19, 20, 24, 25, 26, 27]
        #multiple twofold axis, one high order main axis
        arrgmt = [
            (3, Dict("axis"=>main_axis[2], "relation"=>"∥")),
            (4, Dict())
        ]
    elseif pgidx in [6, 7, 8]
        #multiple twofold axis

        for i in keys(order2_axis)
            j = 1
            while j <= length(order2_axis[i])
                for k = 1:j-1
                    if sort!(order2_axis[i][j][2]) == sort!(order2_axis[i][k][2])
                        deleteat!(order2_axis[i], j)
                        @goto next_j
                    end
                end
                j += 1
                @label next_j
            end
        end

        arrgmt = [
            (3, Dict("axis"=>axis, "relation"=>"∥"))
            for axes in values(order2_axis)
            for (axis, r_iso_type) in axes
        ]
        push!(arrgmt, (4, Dict()))
    elseif pgidx in [10, 11, 17, 22, 23]
        #one high order main axis, improper group
        arrgmt = [
            (2, Dict("axis"=>main_axis[2], "relation"=>"∥")),
            (3, Dict("axis"=>main_axis[2], "relation"=>"∥")),
            (4, Dict())
        ]
    elseif pgidx in [4, 5]
        #one twofold axis, improper group
        arrgmt = [
            (2, Dict("axis"=>main_axis[2], "relation"=>"∥")),
            (3, Dict("axis"=>main_axis[2], "relation"=>"⟂")),
            (3, Dict("axis"=>main_axis[2], "relation"=>"∥")),
            (4, Dict())
        ]
     elseif pgidx == 2
        #-1
        arrgmt = [
            (2, Dict("axis"=>[0, 0, 0])),
            (3, Dict("axis"=>[0, 0, 0])),
            (4, Dict())
        ]
    elseif pgidx in [9, 16, 21]
        #one high order main axis, proper group
        arrgmt = [
            (1, Dict("axis"=>main_axis[2], "relation"=>"∥")),
            (2, Dict("axis"=>main_axis[2], "relation"=>"∥")),
            (3, Dict("axis"=>main_axis[2], "relation"=>"∥")),
            (4, Dict())
        ]
    elseif pgidx == 3
        #one twofold axis, proper group
        arrgmt = [
            (1, Dict("axis"=>main_axis[2], "relation"=>"∥")),
            (2, Dict("axis"=>main_axis[2], "relation"=>"∥")),
            (3, Dict("axis"=>main_axis[2], "relation"=>"⟂")),
            (3, Dict("axis"=>main_axis[2], "relation"=>"∥")),
            (4, Dict())
        ]
    elseif pgidx == 1
        #1
        arrgmt = [
            (1, Dict("axis"=>[0, 0, 0])),
            (2, Dict("axis"=>[0, 0, 0])),
            (3, Dict("axis"=>[0, 0, 0])),
            (4, Dict())
        ]
    end
    return arrgmt
end

end #(module MagneticMoment)

########################################################kSpace########################################################
######################################################################################################################
module kSpace
import ...kSpace as ksp
using LinearAlgebra
using ...MathTool
using ..SpinPointGroupDataBase
using ..SpinPointSymmetry
using ..GroupInformation

export get_symmetry_in_kspace
export get_high_symmetry_klist_of_nontrivial_spin_point_group
export get_high_symmetry_klist_of_full_spin_point_group
export get_high_symmetry_klist_of_spin_point_group
export get_k_invariant_under_group
export get_k_invariant_under_symmetry_operation
export get_symmetry_in_kspace
export get_little_group_of_k
export is_k_invariant_under_group
export is_k_invariant_under_symmetry_operation
export is_k_equivalent_under_group

ktype = Matrix{Float64}
ksymtype = Matrix{Float64}

"""
get_symmetry_in_kspace(rencoded::Union{spg_encodetype, spg_magnonbasis_encodetype}) -> ksymtype
get_symmetry_in_kspace(rsym::spg_symtype) -> ksymtype
get_symmetry_in_kspace(rsym::spg_magnonbasis_symtype) -> ksymtype
"""
function get_symmetry_in_kspace(rencoded::Union{spg_encodetype, spg_magnonbasis_encodetype})::ksymtype
    rsym = spg_decode_symmetry(rencoded)
    return get_symmetry_in_kspace(rsym)
end
function get_symmetry_in_kspace(rsym::spg_symtype)::ksymtype
    ksym = transpose(inv(rsym.r))*det(rsym.s)
    return ksym
end
function get_symmetry_in_kspace(rsym::spg_magnonbasis_symtype)::ksymtype
    ksym = transpose(inv(rsym.r))*rsym.T
    return ksym
end

"""
get_high_symmetry_klist_of_nontrivial_spin_point_group(idx::Int64) -> Dict{Int64, Vector{Any}}
"""
function get_high_symmetry_klist_of_nontrivial_spin_point_group(idx::Int64)::Vector{Matrix{Float64}}
    if !(idx in 1:598) error("Error occurs at get_high_symmetry_klist_of_nontrivial_spin_point_group: nontrivial spin point group index out of range [1, 598]") end
    klist = Vector{Matrix{Float64}}()
    for ks in nspg_high_symmetry_klist[idx]
        push!(klist, ks[1])
    end
    return klist
end


"""
get_high_symmetry_klist_of_spin_point_group(idx::Int64) -> Dict{Int64, Vector{Any}}
"""
function get_high_symmetry_klist_of_spin_point_group(idx::Int64)::Vector{Matrix{Float64}}
    if !(idx in 1:1263) error("Error occurs at get_high_symmetry_klist_of_full_spin_point_group: full spin point group index out of range [1, 1263]") end
    klist = Vector{Matrix{Float64}}()
    for ks in spg_high_symmetry_klist[idx]
        push!(klist, ks[1])
    end
    return klist
end

#=
"""
get_high_symmetry_klist_of_spin_point_group(idx::Int64) -> Dict{Int64, Vector{Any}}
"""
function get_high_symmetry_klist_of_spin_point_group(spgidx::Int64)::Dict{Int64, Vector{Any}}
    if !(spgidx in 1:1263) error("Error occurs at get_high_symmetry_klist_of_spin_point_group: full spin point group index out of range [1, 1263]") end
    group = get_full_symmetry_operations_of_spin_point_group(spgidx, encode=true)
    subgroups = get_full_subgroups_of_spin_point_group(spgidx, encode=true)
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
=#

"""
get_k_invariant_under_nontrivial_spin_point_group(nspgidx::Int64) -> Vector{ktype}
"""
function get_k_invariant_under_nontrivial_spin_point_group(nspgidx::Int64)::Vector{ktype}
    rsyms = get_full_symmetry_operations_of_nontrivial_spin_point_group(nspgidx, encode=false)
    return get_k_invariant_under_group(rsyms)
end

"""
get_k_invariant_under_spin_point_group(spgidx::Int64) -> Vector{ktype}
"""
function get_k_invariant_under_spin_point_group(spgidx::Int64)::Vector{ktype}
    rsyms = get_full_symmetry_operations_of_spin_point_group(spgidx, encode=false)
    return get_k_invariant_under_group(rsyms)
end

"""
get_k_invariant_under_group(rencodeds::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}}) -> Vector{ktype}
get_k_invariant_under_group(rsyms::Vector{<:Union{spg_symtype, spg_magnonbasis_symtype}}) -> Vector{ktype}
"""
function get_k_invariant_under_group(rencodeds::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}})::Vector{ktype}
    rsyms = [spg_decode_symmetry(rencoded) for rencoded in rencodeds]
    return get_k_invariant_under_group(rsyms)
end
function get_k_invariant_under_group(rsyms::Vector{<:Union{spg_symtype, spg_magnonbasis_symtype}})::Vector{ktype}
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
get_k_invariant_under_symmetry_operation(rencoded::spg_encodetype) -> Vector{ktype}
get_k_invariant_under_symmetry_operation(rsym::spg_symtype) -> Vector{ktype}
"""
function get_k_invariant_under_symmetry_operation(rencoded::spg_encodetype)::Vector{ktype}
    rsym = spg_decode_symmetry(rencoded)
    return get_k_invariant_under_symmetry_operation(rsym)
end
function get_k_invariant_under_symmetry_operation(rsym::spg_symtype)::Vector{ktype}
    ksym = get_symmetry_in_kspace(rsym)
    return ksp.get_k_invariant_under_symmetry_operation(ksym)
end

"""
get_little_group_of_k(rencodeds::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}}; k::ktype, encode::Bool) -> Vector{<:Union{spg_symtype, spg_encodetype, spg_magnonbasis_symtype, spg_magnonbasis_encodetype}}
get_little_group_of_k(rsyms::Vector{<:Union{spg_symtype, spg_magnonbasis_symtype}}; k::ktype, encode::Bool) -> Vector{<:Union{spg_symtype, spg_encodetype, spg_magnonbasis_symtype, spg_magnonbasis_encodetype}}
#given point group and k, find little group at k
"""
function get_little_group_of_k(rencodeds::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}}; k::ktype, encode::Bool)::Vector{<:Union{spg_symtype, spg_encodetype, spg_magnonbasis_symtype, spg_magnonbasis_encodetype}}
    rsyms = [spg_decode_symmetry(rencoded) for rencoded in rencodeds]
    return get_little_group_of_k(rsyms, k=k, encode=encode)
end
function get_little_group_of_k(rsyms::Vector{<:Union{spg_symtype, spg_magnonbasis_symtype}}; k::ktype, encode::Bool)::Vector{<:Union{spg_symtype, spg_encodetype, spg_magnonbasis_symtype, spg_magnonbasis_encodetype}}
    little_group = []
    for idx in eachindex(rsyms)
        ksym = get_symmetry_in_kspace(rsyms[idx])
        if ksp.is_k_invariant_under_symmetry_operation(k, ksym) push!(little_group, rsyms[idx]) end
    end
    if encode little_group = [spg_encode_symmetry(rsym) for rsym in little_group] end
    return little_group
end


"""
is_k_invariant_under_group(k::ktype, rencodeds::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}}) -> Bool
is_k_invariant_under_group(k::ktype, rsyms::Vector{<:Union{spg_symtype, spg_magnonbasis_symtype}}) -> Bool
"""
function is_k_invariant_under_group(k::ktype, rencodeds::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}})::Bool
    rsyms = [spg_decode_symmetry(rencoded) for rencoded in rencodeds]
    return is_k_invariant_under_group(k, rsyms)
end
function is_k_invariant_under_group(k::ktype, rsyms::Vector{<:Union{spg_symtype, spg_magnonbasis_symtype}})::Bool
    ksyms = [get_symmetry_in_kspace(rsym) for rsym in rsyms]
    return ksp.is_k_invariant_under_group(k, ksyms)
end


"""
is_k_invariant_under_symmetry_operation(k::ktype, rencoded::Union{spg_encodetype, spg_magnonbasis_encodetype}) -> Bool
is_k_invariant_under_symmetry_operation(k::ktype, rsym::Union{spg_symtype, spg_magnonbasis_symtype}) -> Bool
"""
function is_k_invariant_under_symmetry_operation(k::ktype, rencoded::Union{spg_encodetype, spg_magnonbasis_encodetype})::Bool
    rsym = spg_decode_symmetry(rencoded)
    return is_k_invariant_under_symmetry_operation(k, rsym)
end
function is_k_invariant_under_symmetry_operation(k::ktype, rsym::Union{spg_symtype, spg_magnonbasis_symtype})::Bool
    ksym = get_symmetry_in_kspace(rsym)
    return ksp.is_k_invariant_under_symmetry_operation(k, ksym)
end

"""
is_k_equivalent_under_group(k1::ktype, k2::ktype, rencodeds::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}}) -> Bool
is_k_equivalent_under_group(k1::ktype, k2::ktype, rsyms::Vector{<:Union{spg_symtype, spg_magnonbasis_symtype}}) -> Bool
"""
function is_k_equivalent_under_group(k1::ktype, k2::ktype, rencodeds::Vector{<:Union{spg_encodetype, spg_magnonbasis_encodetype}})::Bool
    rsyms = [spg_decode_symmetry(rencoded) for rencoded in rencodeds]
    return is_k_equivalent_under_group(k1, k2, rsyms)
end
function is_k_equivalent_under_group(k1::ktype, k2::ktype, rsyms::Vector{<:Union{spg_symtype, spg_magnonbasis_symtype}})::Bool
    ksyms = [get_symmetry_in_kspace(rsym) for rsym in rsyms]
    return ksp.is_k_equivalent_under_group(k1, k2, ksyms)
end

end#(module kSpace)

########################################################rSpace########################################################
######################################################################################################################
module rSpace
using LinearAlgebra
using ..SpinPointGroupDataBase
using ..GroupInformation

export generate_magnetic_structure_given_a_site_and_spin_point_group

"""
generate_magnetic_structure_given_a_site_and_nontrivial_spin_point_group(site0, idx::Int64; tor=1e-3)
Input
    site0 = (atom, frac_coord, lattice_magmom)
Return
    site = [(atom, frac_coord, lattice_magmoms, cartesian_magmoms), ...]
    *note lattice vectors and cartesian vectors of magmoms are generally different from realspace lattice and cartesian vectors
"""
function generate_magnetic_structure_given_a_site_and_spin_point_group(site0, idx::Int64; tor=1e-3)
    syms = get_full_symmetry_operations_of_spin_point_group(idx, encode=false)
    sites = [site0]
    for sym in syms
        nsite = (site0[1], mod.(round.(sym.r*site0[2], digits=5),1), sym.s*site0[3])
        for site in sites
            if !(site[1] == nsite[1] && all(abs.(round.(mod.(site[2]-nsite[2], 1), digits=5)) .< tor)) continue end
            if all(abs.(round.(site[3]-nsite[3], digits=5)) .< tor)
                @goto next
            else
                return []
            end
        end
        push!(sites, nsite)
        @label next
    end
    nspgidx = spg_data[idx].nspg
    Bholo = nspg_data[nspgidx].Bholohedry
    if Bholo in [SpinPointGroupDataBase.LaueHolohedry.TRIGO, SpinPointGroupDataBase.LaueHolohedry.HEXA]
        lat = transpose([1.0 0.0 0.0;-0.5 sqrt(3)/2 0.0;0.0 0.0 1.0])
        sites = [(site[1], site[2], site[3], round.(lat*site[3], digits=5)) for site in sites]
    end
    return sites
end

end#(module rSpace)

using .SpinPointGroupDataBase
using .SpinPointSymmetry
using .GroupInformation
using .GroupRepresentation
using .MagneticMoment
using .kSpace
using .rSpace
end #(module SpinPointGroup)
