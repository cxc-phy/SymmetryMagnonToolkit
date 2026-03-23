try
    typeof(MathTool)
catch UndefVarError
    include("../math-tool/mathtool.jl")
end
#include("spacesymmetry.jl")
######################################################################################################################
##################################################SpinSpaceSymmetry###################################################
######################################################################################################################
module SpinSpaceSymmetry
using LinearAlgebra
using ..MathTool
import ..PointSymmetry as pt
import ..SpaceSymmetry as spa
import ..SpinPointSymmetry as spt

export spinspacesym
export magspasym
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
export is_symmetry_unitary
export change_symmetry_from_spin_basis_to_magnon_basis
export isa_group

spinspasym = NamedTuple{(:r, :t, :s), Tuple{Matrix{Int64}, Vector{Float64}, Matrix{Float64}}}
magspasym  = NamedTuple{(:r, :t, :s, :T), Tuple{Matrix{Int64}, Vector{Float64}, Matrix{Float64}, Int64}}

#=
spinspasym = Union{
    NamedTuple{(:r, :t, :s), Tuple{Matrix{Int64}, Vector{Float64}, Matrix{Int64}}},
    NamedTuple{(:r, :t, :s), Tuple{Matrix{Int64}, Vector{Float64}, Matrix{Float64}}}
}
magspasym = Union{
    NamedTuple{(:r, :t, :s, :T), Tuple{Matrix{Int64}, Vector{Float64}, Matrix{Int64}, Int64}},
    NamedTuple{(:r, :t, :s, :T), Tuple{Matrix{Int64}, Vector{Float64}, Matrix{Float64}, Int64}}
}
=#

"""
get_symmetry_in_kspace(rsym::spinspasym)
"""
function get_symmetry_in_kspace(rsym::spinspasym)
    ksym = transpose(inv(rsym.r))*det(rsym.s)
    return ksym
end
function get_symmetry_in_kspace(rsym::magspasym)
    ksym = transpose(inv(rsym.r))*rsym.T
    return ksym
end

"""
get_xyz_form_of_symmetry(sym::spinspasym)
"""
function get_xyz_form_of_symmetry(sym::spinspasym)
    form1 = spa.get_xyz_form_of_symmetry((r=sym.r, t=sym.t))
    form2 = replace(pt.get_xyz_form_of_point_symmetry(Int64.(sym.s)), "x"=>"mx","y"=>"my", "z"=>"mz")
    form = string(form1, form2[end-3:end], "\n", form2[1:end-4])
    return form
end

"""
get_notation_of_symmetry(sym::spinspasym)
"""
function get_notation_of_symmetry(sym::spinspasym)
    R = pt.get_notation_of_point_symmetry(sym.r)
    t = join(sym.t, " ")
    S = pt.get_notation_of_point_symmetry(Int64.(sym.s))
    if R == nothing || S == nothing return nothing end
    notation = string("{", S, "||", R, "|", t, "}")
    return notation
end

"""
is_Identity(syms::spinspasym) -> Bool
"""
function is_Identity(sym::spinspasym)::Bool
    Isym = (r=[1 0 0;0 1 0;0 0 1], t=[0.0, 0.0, 0.0], s=[1.0 0 0;0 1 0;0 0 1])
    return is_symmetry_same(Isym, sym, modtrans=false)
end
function is_Identity(sym::magspasym)::Bool
    Isym = (r=[1 0 0;0 1 0;0 0 1], t=[0.0, 0.0, 0.0], s=[1.0 0;0 1], T=1)
    return is_symmetry_same(Isym, sym, modtrans=false)
end

"""
is_Identity_contained(syms::Vector{spinspasym}) -> Bool
"""
function is_Identity_contained(syms::Vector{spinspasym})::Bool
    Isym = (r=[1 0 0;0 1 0;0 0 1], t=[0.0, 0.0, 0.0], s=[1.0 0 0;0 1 0;0 0 1])
    return is_symmetry_contained(Isym, syms, modtrans=false)
end
function is_Identity_contained(syms::Vector{magspasym})::Bool
    Isym = (r=[1 0 0;0 1 0;0 0 1], t=[0.0, 0.0, 0.0], s=[1.0 0;0 1], T=1)
    return is_symmetry_contained(Isym, syms, modtrans=false)
end

"""
is_symmetry_contained(sym::spinspasym, syms::Vector{spinspasym}; modtrans::Bool) -> Bool
"""
function is_symmetry_contained(sym::T, syms::Vector{T}; modtrans::Bool)::Bool where {T<:Union{magspasym, spinspasym}}
    for ele in syms
        if is_symmetry_same(ele, sym, modtrans=modtrans) return true end
    end
    return false
end

"""
is_symmetry_same(sym1::spinspasym, sym2::spinspasym; modtrans::Bool) -> Bool
"""
function is_symmetry_same(sym1::spinspasym, sym2::spinspasym; modtrans::Bool, tor=1e-3)::Bool
    if modtrans
        relative_t = mod.(round.(sym1.t-sym2.t, digits=3), 1)
    else
        relative_t = round.(sym1.t-sym2.t, digits=3)
    end
    relative_s = round.(sym1.s-sym2.s, digits=3)
    if sym1.r == sym2.r && all(abs.(relative_t) .< tor) && all(abs.(relative_s) .< tor)
        return true
    else
        return false
    end
end
function is_symmetry_same(sym1::magspasym, sym2::magspasym; modtrans::Bool, tor=1e-3)::Bool
    if modtrans
        relative_t = mod.(round.(sym1.t-sym2.t, digits=3), 1)
    else
        relative_t = round.(sym1.t-sym2.t, digits=3)
    end
    relative_s = round.(sym1.s-sym2.s, digits=3)
    if sym1.r == sym2.r && all(abs.(relative_t) .< tor) && all(abs.(relative_s) .< tor) && sym1.T == sym2.T
        return true
    else
        return false
    end
end

"""
get_inverse_of_symmetry(sym::spinspasym; modtrans::Bool) -> spinspasym
"""
function get_inverse_of_symmetry(sym::spinspasym; modtrans::Bool)::spinspasym
    r_inv = get_inverse_of_matrix(sym.r)
    s_inv = get_inverse_of_matrix(sym.s)
    if modtrans
        sym_inv = (r=r_inv, t=mod.(-r_inv*sym.t, 1), s=s_inv)
    else
        sym_inv = (r=r_inv, t=-r_inv*sym.t, s=s_inv)
    end
    return sym_inv
end
function get_inverse_of_symmetry(sym::magspasym; modtrans::Bool)::magspasym
    r_inv = get_inverse_of_matrix(sym.r)
    s_inv = get_inverse_of_matrix(sym.s)
    if modtrans
        sym_inv = (r=r_inv, t=mod.(-r_inv*sym.t, 1), s=s_inv, T=sym.T)
    else
        sym_inv = (r=r_inv, t=-r_inv*sym.t, s=s_inv, T=sym.T)
    end
    return sym_inv
end

"""
composite_symmetries(syms::spinspasym; order::Integer, modtrans::Bool) -> spinspasym
composite_symmetries(syms::spinspasym...; modtrans::Bool) -> spinspasym
"""
function composite_symmetries(sym::spinspasym; order::Integer=1, modtrans::Bool)::spinspasym
    if order == 1
        return sym
    elseif order > 1
        return composite_symmetries([sym for i = 1:order]..., modtrans=modtrans)
    else
        error("Error occurs at composite_symmetries: order should be a positive integer")
    end
end
function composite_symmetries(syms::spinspasym...; modtrans::Bool)::spinspasym
    if length(syms) == 1
        return syms[1]
    elseif length(syms) == 2
        if modtrans
            return (r=syms[1].r*syms[2].r, t=mod.(syms[1].t+syms[1].r*syms[2].t,1), s=syms[1].s*syms[2].s)
        else
            return (r=syms[1].r*syms[2].r, t=syms[1].t+syms[1].r*syms[2].t, s=syms[1].s*syms[2].s)
        end
    else
        return composite_symmetries(composite_symmetries(syms[1], syms[2], modtrans = modtrans), syms[3:end]..., modtrans = modtrans)
    end
end
function composite_symmetries(sym::magspasym; order::Integer=1, modtrans::Bool)::magspasym
    if order == 1
        return sym
    elseif order > 1
        return composite_symmetries([sym for i = 1:order]..., modtrans=modtrans)
    else
        error("Error occurs at composite_symmetries: order should be a positive integer")
    end
end
function composite_symmetries(syms::magspasym...; modtrans::Bool)::magspasym
    if length(syms) == 1
        return syms[1]
    elseif length(syms) == 2
        if modtrans
            nsym = (r=syms[1].r*syms[2].r, t=mod.(syms[1].t+syms[1].r*syms[2].t,1), s=syms[1].s*syms[2].s, T=syms[1].T*syms[2].T)
        else
            nsym = (r=syms[1].r*syms[2].r, t=syms[1].t+syms[1].r*syms[2].t, s=syms[1].s*syms[2].s, T=syms[1].T*syms[2].T)
        end
        return nsym
    else
        return composite_symmetries(composite_symmetries(syms[1], syms[2], modtrans=modtrans), syms[3:end]..., modtrans=modtrans)
    end
end

"""
conjugate_symmetry(sym::magspasym, conjg_sym::magspasym; modtrans::Bool)->magspasym
"""
function conjugate_symmetry(sym::magspasym, conjg_sym::magspasym; modtrans::Bool)::magspasym
    return composite_symmetries(get_inverse_of_symmetry(conjg_sym, modtrans=modtrans), sym, conjg_sym, modtrans=modtrans)
end

"""
is_symmetry_unitary(sym::spinspasym) -> Bool
is_symmetry_unitary(sym::magspasym) -> Bool
"""
function is_symmetry_unitary(sym::spinspasym)::Bool
    if det(sym.s) > 0
        return true
    else
        return false
    end
end
function is_symmetry_unitary(sym::magspasym)::Bool
    if sym.T == 1
        return true
    else
        return false
    end
end

"""
change_symmetry_from_spin_basis_to_magnon_basis(sym::spinspasym; SO2::Vector{<:Real}) -> magspasym
"""
function change_symmetry_from_spin_basis_to_magnon_basis(sym::spinspasym; SO2::Vector{<:Real})::magspasym
    nsym = spt.change_symmetry_from_spin_basis_to_magnon_basis((r=sym.r, s=sym.s), SO2=SO2)
    nsym = (r=nsym.r, t=sym.t, s=real.(nsym.s),T=nsym.T)
    return nsym
end

"""
isa_group(syms::Vector{T}) -> Bool where {T<:Union{spinspasym, magspasym}}
check if a set forms a group
"""
function isa_group(syms::Vector{T})::Bool where {T<:Union{spinspasym, magspasym}}
    if !is_Identity_contained(syms) return false end
    for ele1 in syms
        ele1_inv = get_inverse_of_symmetry(ele1, modtrans=false)
        if !is_symmetry_contained(ele1_inv, syms, modtrans=true)
            return false
        end
        for ele2 in syms
            ele3 = composite_symmetries(ele1, ele2, modtrans=false)
            if !is_symmetry_contained(ele3, syms, modtrans=true)
                return false
            end
        end
    end
    return true
end

"""
get_coset_representatives_wrt_translation_group(syms::Vector{sg.sg_symtype}) -> Vector{sg.sg_symtype}
"""
function get_coset_representatives_wrt_translation_group(syms::Vector{spinspasym}; tor=1e-5)::Vector{spinspasym}
    cosetrep = spinspasym[]
    for sym in syms
        for rep in cosetrep
            if sym.r == rep.r && all(abs.(round.(sym.s-rep.s, digits=5)) .< tor) @goto next_sym end
        end
        push!(cosetrep, sym)
        @label next_sym
    end
    return cosetrep
end
function get_coset_representatives_wrt_translation_group(syms::Vector{magspasym}; tor=1e-5)::Vector{magspasym}
    cosetrep = magspasym[]
    for sym in syms
        for rep in cosetrep
            if sym.r == rep.r && all(abs.(round.(sym.s-rep.s, digits=5)) .< tor) && sym.T == rep.T @goto next_sym end
        end
        push!(cosetrep, sym)
        @label next_sym
    end
    return cosetrep
end


end#(module SpinSpaceSymmetry)
