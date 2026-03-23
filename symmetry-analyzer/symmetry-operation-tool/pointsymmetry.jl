try
    typeof(MathTool)
catch UndefVarError
    include("../math-tool/mathtool.jl")
end
######################################################################################################################
####################################################PointSymmetry#####################################################
######################################################################################################################
module PointSymmetry
using LinearAlgebra
using ..MathTool

export ptsym
export get_symmetry_in_kspace
export get_xyz_form_of_point_symmetry
export get_notation_of_point_symmetry
export is_Identity
export is_Identity_contained
export is_symmetry_contained
export is_symmetry_same
export get_inverse_of_symmetry
export composite_symmetries
export isa_group


export get_symmetry_order_of_point_symmetry
export get_isometry_type_of_point_symmetry
export get_isometry_axis_of_point_symmetry
export get_proper_part_of_rotation_of_point_symmetry
export get_sense_of_rotation_of_point_symmetry

ptsym = Matrix{Int64}

"""
get_symmetry_in_kspace(rsym::ptsym)
"""
function get_symmetry_in_kspace(rsym::ptsym)
    ksym = transpose(inv(rsym))
    return ksym
end

"""
get_xyz_form_of_point_symmetry(sym::ptsym; tor=1e-3)
"""
function get_xyz_form_of_point_symmetry(sym::ptsym; tor=1e-3)
    form = ["", "", "", ""]
    xyz = ["x", "y", "z"]
    for row in 1:3
        for idx in 1:3
            if sym[row, idx] == 1
                if length(strip(form[row])) != 0 form[row] = string(form[row], "+") end
                form[row] = string(form[row], xyz[idx])
            elseif sym[row, idx] == -1
                form[row] = string(form[row], "-")
                form[row] = string(form[row], xyz[idx])
            end

        end
    end
    if det(sym) > tor
        form[4] = string(form[4], "+1")
    else
        form[4] = string(form[4], "-1")
    end
    form = join(form, ", ")
    return form
end

"""
get_notation_of_point_symmetry(sym::ptsym)
operation    -6 -4 -3 -2 -1  1  2  3  4  6
index         1  2  3  4  5  6  7  8  9 10
"""
function get_notation_of_point_symmetry(sym::ptsym)
    iso_type = get_isometry_type_of_point_symmetry(sym)
    iso_axis = get_isometry_axis_of_point_symmetry(sym)
    rot_sense = get_sense_of_rotation_of_point_symmetry(sym)
    notation = nothing
    if iso_type == 1
        if rot_sense == 1
            notation = string("C_", iso_axis, "(π/3)", "I")
        elseif rot_sense == -1
            notation = string("C_", iso_axis, "(5π/3)", "I")
        end
    elseif iso_type == 2
        if rot_sense == 1
            notation = string("C_", iso_axis, "(π/2)", "I")
        elseif rot_sense == -1
            notation = string("C_", iso_axis, "(3π/2)", "I")
        end
    elseif iso_type == 3
        if rot_sense == 1
            notation = string("C_", iso_axis, "(2π/3)", "I")
        elseif rot_sense == -1
            notation = string("C_", iso_axis, "(4π/3)", "I")
        end
    elseif iso_type == 4
        notation = string("C_", iso_axis, "(π)", "I")
    elseif iso_type == 5
        notation = "I"
    elseif iso_type == 6
        notation = "E"
    elseif iso_type == 7
        notation = string("C_", iso_axis, "(π)")
    elseif iso_type == 8
        if rot_sense == 1
            notation = string("C_", iso_axis, "(2π/3)")
        elseif rot_sense == -1
            notation = string("C_", iso_axis, "(4π/3)")
        end
    elseif iso_type == 9
        if rot_sense == 1
            notation = string("C_", iso_axis, "(π/2)")
        elseif rot_sense == -1
            notation = string("C_", iso_axis, "(3π/2)")
        end
    elseif iso_type == 10
        if rot_sense == 1
            notation = string("C_", iso_axis, "(π/3)")
        elseif rot_sense == -1
            notation = string("C_", iso_axis, "(5π/3)")
        end
    end
    return notation
end

"""
is_Identity(sym::ptsym) -> Bool
"""
function is_Identity(sym::ptsym)::Bool
    Isym = [1 0 0;0 1 0;0 0 1]
    return is_symmetry_same(Isym, sym)
end

"""
is_Identity_contained(syms::Vector{ptsym}) -> Bool
"""
function is_Identity_contained(syms::Vector{ptsym})::Bool
    Isym = [1 0 0;0 1 0;0 0 1]
    return is_symmetry_contained(Isym, syms)
end

"""
is_symmetry_contained(sym::ptsym, syms::Vector{ptsym}) -> Bool
"""
function is_symmetry_contained(sym::ptsym, syms::Vector{ptsym})::Bool
    for ele in syms
        if is_symmetry_same(ele, sym) return true end
    end
    return false
end

"""
is_symmetry_same(sym1::ptsym, sym2::ptsym) -> Bool
"""
function is_symmetry_same(sym1::ptsym, sym2::ptsym)::Bool
    if sym1 == sym2
        return true
    else
        return false
    end
end

"""
get_inverse_of_symmetry(sym::ptsym) -> ptsym
"""
function get_inverse_of_symmetry(sym::ptsym)::ptsym
    sym_inv = get_inverse_of_matrix(sym)
    return sym_inv
end

"""
composite_symmetries(sym::ptsym; order::Integer) -> ptsym
composite_symmetries(syms::ptsym...) -> ptsym
"""
function composite_symmetries(sym::ptsym; order::Integer)::ptsym
    if order == 1
        return sym
    elseif order > 1
        return composite_symmetries([sym for i = 1:order]...)
    else
        error("Error occurs at composite_symmetries: order should be a positive integer")
    end
end
function composite_symmetries(syms::Vector{ptsym})::ptsym
    return composite_symmetries(syms...)
end
function composite_symmetries(syms::ptsym...)::ptsym
    if length(syms) == 1
        return syms[1]
    elseif length(syms) == 2
        return syms[1]*syms[2]
    else
        return composite_symmetries(composite_symmetries(syms[1], syms[2]), syms[3:end]...)
    end
end

"""
isa_group(syms::Vector{ptsym}) -> Bool
check if a set forms a group
"""
function isa_group(syms::Vector{ptsym})::Bool
    if !is_Identity_contained(syms) return false end
    for ele1 in syms
        ele1_inv = get_inverse_of_symmetry(ele1)
        if !is_symmetry_contained(ele1_inv, syms)
            return false
        end
        for ele2 in syms
            ele3 = composite_symmetries(ele1, ele2)
            if !is_symmetry_contained(ele3, syms)
                return false
            end
        end
    end
    return true
end

"""
get_symmetry_order_of_point_symmetry(sym::ptsym) -> Union{Int64, Nothing}
#get order n of a point symmetry operation, sym^n=I
"""
function get_symmetry_order_of_point_symmetry(sym::ptsym)::Union{Int64, Nothing}

    iso_type = get_isometry_type_of_point_symmetry(sym)
    order = nothing
    if iso_type == 1 || iso_type == 3 || iso_type == 10
        order = 6
    elseif iso_type == 2 || iso_type == 9
        order = 4
    elseif iso_type == 4 || iso_type == 5 || iso_type == 7
        order = 2
    elseif iso_type == 6
        order = 1
    elseif iso_type == 8
        order = 3
    end

    #if order == nothing println("Warning: symmetry order of point symmetry not found") end
    return order
end

"""
get_isometry_type_of_point_symmetry(sym::ptsym) -> Union{Int64, Nothing}
isometry_type
operation    -6 -4 -3 -2 -1  1  2  3  4  6
order         6  4  6  2  2  1  2  3  4  6
trace        -2 -1  0  1 -3  3 -1  0  1  2
determinant  -1 -1 -1 -1 -1  1  1  1  1  1
index         1  2  3  4  5  6  7  8  9 10


for sym belonging to 32 point groups
transpose(sym)*sym in [
    [1  0 0;  0 1 0; 0 0 1],
    [1 -1 0; -1 2 0; 0 0 1],
    [2 -1 0; -1 1 0; 0 0 1]
]
"""
function get_isometry_type_of_point_symmetry(sym::ptsym)::Union{Int64, Nothing}
    iso_type = nothing
    if det(sym) == -1
        if tr(sym) == -2
            iso_type = 1
            order = 6
        elseif tr(sym) == -1
            iso_type = 2
            order = 4
        elseif tr(sym) == 0
            iso_type = 3
            order = 6
        elseif tr(sym) == 1
            iso_type = 4
            order = 2
        elseif tr(sym) == -3
            iso_type = 5
            order = 2
        end
    elseif det(sym)== 1
        if tr(sym) == 3
            iso_type = 6
            order = 1
        elseif tr(sym) == -1
            iso_type = 7
            order = 2
        elseif tr(sym) == 0
            iso_type = 8
            order = 3
        elseif tr(sym) == 1
            iso_type = 9
            order = 4
        elseif tr(sym) == 2
            iso_type = 10
            order = 6
        end
    end
    return iso_type
end

"""
get_isometry_axis_of_point_symmetry(sym::ptsym) -> Union{Vector{Int64}, Nothing}
#formula: Y(W,k)=W^1+W^2+...+W^k
#method
    if sym = rotation, Y(sym, order)
    if sym = roto-inversion / reflection, Y(-sym, order)
"""
function get_isometry_axis_of_point_symmetry(sym::ptsym; tor=1e-3)::Union{Vector{Int64}, Nothing}
    proper = get_proper_part_of_point_symmetry(sym)
    order = get_symmetry_order_of_point_symmetry(proper)
    axis = nothing
    if order == 1
        axis = [0, 0, 0]
    elseif order in [2, 3, 4, 6]
        Y = sum(proper^i for i = 1:order)
        sum_column = sum(Y[:, i]*[sign(ele) for ele in Y[:,i] if sign(ele)!=0][1] for i = 1:3 if any(Y[:,i].!=0))
        axis = sum_column/gcd(sum_column)
    end

    #=
    if axis == nothing
        println("Warning: isometry axis of point symmetry not found")
    end
    =#
    return axis
end

"""
get_proper_part_of_point_symmetry(sym::ptsym) -> ptsym
#get proper part of a pointsymmetry operation
"""
function get_proper_part_of_point_symmetry(sym::ptsym)::ptsym
    proper = Int64.(det(sym)*sym)
    return proper
end

"""
get_sense_of_rotation_of_point_symmetry(sym::ptsym) -> Union{Int64, Nothing}
#formula: det(z),z=[u|x|det(sym)*(sym*x)]
"""
function get_sense_of_rotation_of_point_symmetry(sym::ptsym)::Union{Int64, Nothing}
    rot_sense = nothing
    iso_type = get_isometry_type_of_point_symmetry(sym)
    if iso_type in [4, 5, 6, 7]
        rot_sense = 0
    elseif iso_type in [1, 2, 3, 8, 9, 10]
        u = get_isometry_axis_of_point_symmetry(sym)
        for j in [1, 2, 3]
            x = sym[:, j]
            z = [u;; x;; det(sym)*(sym*x)]
            rot_sense = Int64.(round(det(z), digits=5))
            if rot_sense in [1, -1]
                break
            else
                rot_sense = nothing
            end
        end
    end

    #if rot_sense == nothing
    #    println("Warning: sense of rotation of point symmetry not found")
    #end
    return rot_sense
end

end#(module PointSymmetry)
