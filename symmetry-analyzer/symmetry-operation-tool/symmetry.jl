include("./pointsymmetry.jl")
include("./spacesymmetry.jl")
include("./spinpointsymmetry.jl")
include("./spinspacesymmetry.jl")
try
    typeof(MathTool)
catch UndefVarError
    include("../math-tool/mathtool.jl")
end
######################################################################################################################
######################################################Symmetry########################################################
######################################################################################################################
module Symmetry
using LinearAlgebra
using ..MathTool
import ..PointSymmetry as pt
import ..SpaceSymmetry as spa
import ..SpinPointSymmetry as spinpt
import ..SpinSpaceSymmetry as spinspa

using .pt:get_symmetry_order_of_point_symmetry
using .pt:get_isometry_axis_of_point_symmetry
using .pt:get_notation_of_point_symmetry
export get_symmetry_order_of_point_symmetry
export get_isometry_axis_of_point_symmetry
export get_notation_of_point_symmetry

export ptsym
export spasym
export spinptsym
export magptsym
export spinspasym
export magspasym

#export get_kxkykz_form_of_symmetry
export get_symmetry_in_kspace
export get_xyz_form_of_symmetry
export get_notation_of_symmetry
export is_Identity
export is_Identity_contained
export is_symmetry_contained
export is_symmetry_same
export get_inverse_of_symmetry
export composite_symmetries
export conjugate_symmetry
export isa_pointsymmetry
export is_symmetry_unitary
export change_symmetry_from_spin_basis_to_magnon_basis
export get_multiplication_table_of_group
export isa_group
export get_coset_representatives_wrt_translation_group


ptsym = pt.ptsym
spasym = spa.spasym
spinptsym = spinpt.spinptsym
magptsym = spinpt.magptsym
spinspasym = spinspa.spinspasym
magspasym  = spinspa.magspasym
Ttype = Union{ptsym, spasym, spinptsym, spinspasym, magptsym, magspasym}
#=
"""
get_kxkykz_form_of_symmetry(rsym::Union{ptsym, spasym, spinptsym, spinspasym})
"""
function get_kxkykz_form_of_symmetry(rsym::Union{ptsym, spasym, spinptsym, spinspasym})
    ksym = get_symmetry_in_kspace(rsym)
    form = get_xyz_form_of_symmetry(ksym)
    form = replace(form, "x"=>"kx", "y"=>"ky", "z"=>"kz")
    return form
end
=#

"""
get_symmetry_in_kspace(rsym::Ttype)
"""
function get_symmetry_in_kspace(rsym::Ttype)
    if isa(rsym, ptsym)
        return pt.get_symmetry_in_kspace(rsym)
    elseif isa(rsym, spasym)
        return spa.get_symmetry_in_kspace(rsym)
    elseif isa(rsym, Union{spinptsym, magptsym})
        return spinpt.get_symmetry_in_kspace(rsym)
    elseif isa(rsym, Union{spinspasym, magspasym})
        return spinspa.get_symmetry_in_kspace(rsym)
    end
end

"""
get_xyz_form_of_symmetry(sym::Ttype)
"""
function get_xyz_form_of_symmetry(rsym::Ttype)
    if isa(rsym, ptsym)
        return pt.get_xyz_form_of_point_symmetry(rsym)[1:end-4]
    elseif isa(rsym, spasym)
        return spa.get_xyz_form_of_symmetry(rsym)
    elseif isa(rsym, Union{spinptsym, magptsym})
        return spinpt.get_xyz_form_of_symmetry(rsym)
    elseif isa(rsym, Union{spinspasym, magspasym})
        return spinspa.get_xyz_form_of_symmetry(rsym)
    end
end

"""
get_notation_of_symmetry(sym::Ttype)
"""
function get_notation_of_symmetry(sym::Ttype)
    if isa(sym, ptsym)
        return pt.get_notation_of_point_symmetry(sym)
    elseif isa(sym, spasym)
        return spa.get_notation_of_symmetry(sym)
    elseif isa(sym, Union{spinptsym, magptsym})
        return spinpt.get_notation_of_symmetry(sym)
    elseif isa(sym, Union{spinspasym, magspasym})
        return spinspa.get_notation_of_symmetry(sym)
    end
end

"""
is_Identity(sym::T) -> Bool where {T<:Ttype}
"""
function is_Identity(sym::Ttype)::Bool
    if isa(sym, ptsym)
        return pt.is_Identity(sym)
    elseif isa(sym, spasym)
        return spa.is_Identity(sym)
    elseif isa(sym, Union{spinptsym, magptsym})
        return spinpt.is_Identity(sym)
    elseif isa(sym, Union{spinspasym, magspasym})
        return spinspa.is_Identity(sym)
    end
end

"""
is_Identity_contained(syms::Vector{T}) -> Bool where {T<:Ttype}
"""
function is_Identity_contained(syms::Vector{T})::Bool where {T<:Ttype}
    if isa(syms[1], ptsym)
        return pt.is_Identity_contained(syms)
    elseif isa(syms[1], spasym)
        return spa.is_Identity_contained(syms)
    elseif isa(syms[1], Union{spinptsym, magptsym})
        return spinpt.is_Identity_contained(syms)
    elseif isa(syms[1], Union{spinspasym, magspasym})
        return spinspa.is_Identity_contained(syms)
    end
end

"""
is_symmetry_contained(sym::T, syms::Vector{T}; modtrans=nothing) -> Bool where {T<:Ttype}
"""
function is_symmetry_contained(sym::T, syms::Vector{T}; modtrans=nothing)::Bool where {T<:Ttype}
    if isa(sym, ptsym)
        return pt.is_symmetry_contained(sym, syms)
    elseif isa(sym, spasym)
        return spa.is_symmetry_contained(sym, syms, modtrans=modtrans)
    elseif isa(sym, Union{spinptsym, magptsym})
        return spinpt.is_symmetry_contained(sym, syms)
    elseif isa(sym, Union{spinspasym, magspasym})
        return spinspa.is_symmetry_contained(sym, syms, modtrans=modtrans)
    end
end

"""
is_symmetry_contained(sym1::T, sym2::T) -> Bool where {T<:Ttype}
"""
function is_symmetry_same(sym1::T, sym2::T; modtrans=nothing)::Bool where {T<:Ttype}
    if isa(sym1, ptsym)
        return pt.is_symmetry_same(sym1, sym2)
    elseif isa(sym1, spasym)
        return spa.is_symmetry_same(sym1, sym2, modtrans=modtrans)
    elseif isa(sym1, Union{spinptsym, magptsym})
        return spinpt.is_symmetry_same(sym1, sym2)
    elseif isa(sym1, Union{spinspasym, magspasym})
        return spinspa.is_symmetry_same(sym1, sym2, modtrans=modtrans)
    end
end

"""
get_inverse_of_symmetry(sym::T; modtrans=nothing) -> T where {T<:Ttype}
"""
function get_inverse_of_symmetry(sym::T; modtrans=nothing)::T where {T<:Ttype}
    if isa(sym, ptsym)
        return pt.get_inverse_of_symmetry(sym)
    elseif isa(sym, spasym)
        return spa.get_inverse_of_symmetry(sym, modtrans=modtrans)
    elseif isa(sym, Union{spinptsym, magptsym})
        return spinpt.get_inverse_of_symmetry(sym)
    elseif isa(sym, Union{spinspasym, magspasym})
        return spinspa.get_inverse_of_symmetry(sym, modtrans=modtrans)
    end
end

"""
composite_symmetries(sym::T; order::Integer) -> T where {T<:Ttype}
composite_symmetries(syms::Vector{T}) -> T where {T<:Ttype}
composite_symmetries(syms::T...) -> T where {T<:Ttype}
"""
function composite_symmetries(sym::Ttype; order::Integer=1, modtrans=nothing)::Ttype# where {T<:Ttype}
    if isa(sym, ptsym)
        return pt.composite_symmetries(sym, order=order)
    elseif isa(sym, spasym)
        return spa.composite_symmetries(sym, order, modtrans=modtrans)
    elseif isa(sym, Union{spinptsym, magptsym})
        return spinpt.composite_symmetries(sym, order=order)
    elseif isa(sym, Union{spinspasym, magspasym})
        return spinspa.composite_symmetries(sym, order=order, modtrans=modtrans)
    end
end
function composite_symmetries(syms::Vector{<:Ttype}; modtrans=nothing)::Ttype #where {T<:Ttype}
    return composite_symmetries(syms..., modtrans=modtrans)
end
function composite_symmetries(syms::Ttype...; modtrans=nothing)::Ttype# where {T<:Ttype}
    if isa(syms[1], ptsym)
        return pt.composite_symmetries(syms...)
    elseif isa(syms[1], spasym)
        return spa.composite_symmetries(syms..., modtrans=modtrans)
    elseif isa(syms[1], Union{spinptsym, magptsym})
        return spinpt.composite_symmetries(syms...)
    elseif isa(syms[1], Union{spinspasym, magspasym})
        return spinspa.composite_symmetries(syms..., modtrans=modtrans)
    end
end

"""
conjugate_symmetry(sym::T, conjg_sym::T; modtrans=nothing) -> T where {T<:Ttype}
"""
function conjugate_symmetry(sym::T, conjg_sym::T; modtrans=nothing)::T where {T<:Ttype}
    if isa(sym, ptsym)
        error("not applied")
        #return pt.conjugate_symmetry(sym, conjg_sym)
    elseif isa(sym, spasym)
        error("not applied")
        #return spa.conjugate_symmetry(sym, conjg_sym, modtrans=modtrans)
    elseif isa(sym, Union{spinptsym, magptsym})
        error("not applied")
        #return spinpt.conjugate_symmetry(sym, conjg_sym)
    elseif isa(sym, Union{spinspasym, magspasym})
        return spinspa.conjugate_symmetry(sym, conjg_sym, modtrans=modtrans)
    end
end

"""
isa_pointsymmetry(sym::spasym) -> Bool
"""
function isa_pointsymmetry(sym::spasym)::Bool
    return spa.isa_pointsymmetry(sym)
end

"""
is_symmetry_unitary(sym::Union{spinptsym, magptsym, spinspasym, magspasym}) -> Bool
"""
function is_symmetry_unitary(sym::Union{spinptsym, magptsym, spinspasym, magspasym})::Bool
    if isa(sym, Union{spinptsym, magptsym})
        return spinpt.is_symmetry_unitary(sym)
    elseif isa(sym, Union{spinspasym, magspasym})
        return spinspa.is_symmetry_unitary(sym)
    end
    return spa.isa_pointsymmetry(sym)
end

"""
change_symmetry_from_spin_basis_to_magnon_basis(sym::Union{spinptsym, spinspasym}) -> Union{magptsym, magspasym}
"""
function change_symmetry_from_spin_basis_to_magnon_basis(sym::Union{spinptsym, spinspasym}; SO2::Vector{<:Real})::Union{magptsym, magspasym}
    if isa(sym, spinptsym)
        return spinpt.change_symmetry_from_spin_basis_to_magnon_basis(sym, SO2=SO2)
    elseif isa(sym, spinspasym)
        return spinspa.change_symmetry_from_spin_basis_to_magnon_basis(sym, SO2=SO2)
    end
end


"""
get_multiplication_table_of_group(syms::Vector{T}) -> Matrix{Int64} where {T<:Ttype}
"""
function get_multiplication_table_of_group(syms::Vector{T})::Matrix{Int64} where {T<:Ttype}
    if isa(syms[1], ptsym)
        return pt.get_multiplication_table_of_group(syms)
    elseif isa(syms[1], spasym)
        return spa.get_multiplication_table_of_group(syms)
    elseif isa(syms[1], Union{spinptsym, magptsym})
        return spinpt.get_multiplication_table_of_group(syms)
    elseif isa(syms[1], Union{spinspasym, magspasym})
        return spinspa.get_multiplication_table_of_group(syms)
    end
end

"""
isa_group(syms::Vector{T}) -> Bool where {T<:Union{ptsym, spasym, spinptsym, spinspasym}}
check if a set forms a group
"""
function isa_group(syms::Vector{T})::Bool where {T<:Ttype}
    if isa(syms[1], ptsym)
        return pt.isa_group(syms)
    elseif isa(syms[1], spasym)
        return spa.isa_group(syms)
    elseif isa(syms[1], Union{spinptsym, magptsym})
        return spinpt.isa_group(syms)
    elseif isa(syms[1], Union{spinspasym, magspasym})
        return spinspa.isa_group(syms)
    end
end

"""
isa_group(syms::Vector{T}) -> Bool where {T<:Union{ptsym, spasym, spinptsym, spinspasym}}
check if a set forms a group
"""
function get_coset_representatives_wrt_translation_group(syms::Vector{T})::Vector{T} where {T<:Union{spasym, spinspasym, magspasym}}
    if isa(syms[1], spasym)
        error("not applied")
    elseif isa(syms[1], Union{spinspasym, magspasym})
        return spinspa.get_coset_representatives_wrt_translation_group(syms)
    end
end

end#(module Symmetry)
