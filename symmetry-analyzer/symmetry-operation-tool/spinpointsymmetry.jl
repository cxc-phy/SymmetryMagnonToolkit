try
    typeof(MathTool)
catch UndefVarError
    include("../math-tool/mathtool.jl")
end
######################################################################################################################
##################################################SpinPointSymmetry###################################################
######################################################################################################################
module SpinPointSymmetry
using LinearAlgebra
using ..MathTool
import ..PointSymmetry as pt


export spinptsym
export magptsym
export get_symmetry_in_kspace
export get_xyz_form_of_symmetry
export get_notation_of_symmetry
export is_Identity
export is_Identity_contained
export is_symmetry_contained
export is_symmetry_same
export get_inverse_of_symmetry
export composite_symmetries
export is_symmetry_unitary
export change_symmetry_from_spin_basis_to_magnon_basis
export get_multiplication_table_of_group
export isa_group

spinptsym = NamedTuple{(:r, :s), Tuple{Matrix{Int64}, Matrix{Float64}}}
magptsym = NamedTuple{(:r, :s, :T), Tuple{Matrix{Int64}, Matrix{Float64}, Int64}}

#=
spinptsym = Union{
    NamedTuple{(:r, :s), Tuple{Matrix{Int64}, Matrix{Int64}}},
    NamedTuple{(:r, :s), Tuple{Matrix{Int64}, Matrix{Float64}}}
}
magptsym = Union{
    NamedTuple{(:r, :s, :T), Tuple{Matrix{Int64}, Matrix{Int64}, Int64}},
    NamedTuple{(:r, :s, :T), Tuple{Matrix{Int64}, Matrix{Float64}, Int64}}
}
=#
"""
get_symmetry_in_kspace(rsym::spinptsym)
"""
function get_symmetry_in_kspace(rsym::spinptsym)
    ksym = transpose(inv(rsym.r))*det(rsym.s)
    return ksym
end

"""
get_xyz_form_of_symmetry(sym::spinptsym)
"""
function get_xyz_form_of_symmetry(sym::spinptsym)
    form1 = pt.get_xyz_form_of_point_symmetry(sym.r)[1:end-4]
    form2 = replace(pt.get_xyz_form_of_point_symmetry(Int64.(sym.s)), "x"=>"mx","y"=>"my", "z"=>"mz")
    form = string(form1, form2[end-3:end], "\n", form2[1:end-4])
    return form
end

"""
get_notation_of_symmetry(sym::spinptsym)
"""
function get_notation_of_symmetry(sym::spinptsym)
    R = pt.get_notation_of_point_symmetry(sym.r)
    t = "0 0 0"
    S = pt.get_notation_of_point_symmetry(Int64.(sym.s))
    if R == nothing || S == nothing return nothing end
    notation = string("{", S, "||", R, "|", t, "}")
    return notation
end

"""
is_Identity(sym::spinptsym) -> Bool
"""
function is_Identity(sym::spinptsym)::Bool
    Isym = (r=[1 0 0;0 1 0;0 0 1], s=[1.0 0 0;0 1 0;0 0 1])
    return is_symmetry_same(Isym, sym)
end
function is_Identity(sym::magptsym)::Bool
    Isym = (r=[1 0 0;0 1 0;0 0 1], s=[1.0 0;0 1], T=1)
    return is_symmetry_same(Isym, sym)
end

"""
is_Identity_contained(syms::Vector{spinptsym}) -> Bool
"""
function is_Identity_contained(syms::Vector{spinptsym})::Bool
    Isym = (r=[1 0 0;0 1 0;0 0 1], s=[1.0 0 0;0 1 0;0 0 1])
    return is_symmetry_contained(Isym, syms)
end
function is_Identity_contained(syms::Vector{magptsym})::Bool
    Isym = (r=[1 0 0;0 1 0;0 0 1], s=[1.0 0;0 1], T=1)
    return is_symmetry_contained(Isym, syms)
end

"""
is_symmetry_contained(sym::spinptsym, syms::Vector{spinptsym}) -> Bool
"""
function is_symmetry_contained(sym::T, syms::Vector{T})::Bool where {T<:Union{magptsym, spinptsym}}
    for ele in syms
        if is_symmetry_same(ele, sym) return true end
    end
    return false
end

"""
is_symmetry_same(sym1::spinptsym, sym2::spinptsym) -> Bool
"""
function is_symmetry_same(sym1::spinptsym, sym2::spinptsym; tor=1e-5)::Bool
    if sym1.r == sym2.r && all(abs.(round.(sym1.s - sym2.s, digits=5)) .< tor)
        return true
    else
        return false
    end
end
function is_symmetry_same(sym1::magptsym, sym2::magptsym; tor=1e-5)::Bool
    if sym1.r == sym2.r && all(abs.(round.(sym1.s - sym2.s, digits=5)) .< tor) && sym1.T == sym2.T
        return true
    else
        return false
    end
end

"""
get_inverse_of_symmetry(sym::spinptsym) -> spinptsym
"""
function get_inverse_of_symmetry(sym::spinptsym)::spinptsym
    r_inv = get_inverse_of_matrix(sym.r)
    s_inv = get_inverse_of_matrix(sym.s)
    sym_inv = (r=r_inv, s=s_inv)
    return sym_inv
end
function get_inverse_of_symmetry(sym::magptsym)::magptsym
    r_inv = get_inverse_of_matrix(sym.r)
    s_inv = get_inverse_of_matrix(sym.s)
    sym_inv = (r=r_inv, s=s_inv, T=sym.T)
    return sym_inv
end

"""
composite_symmetries(sym::spinptsym; order::Integer) -> spinptsym
composite_symmetries(syms::spinptsym...) -> spinptsym
"""
function composite_symmetries(sym::spinptsym; order::Integer)::spinptsym
    if order == 1
        return sym
    elseif order > 1
        return composite_symmetries([sym for i = 1:order]...)
    else
        error("Error occurs at composite_symmetries: order should be a positive integer")
    end
end
function composite_symmetries(syms::spinptsym...)::spinptsym
    if length(syms) == 1
        return syms[1]
    elseif length(syms) == 2
        nsym = (r=syms[1].r*syms[2].r, s=syms[1].s*syms[2].s)
        return nsym
    else
        return composite_symmetries(composite_symmetries(syms[1], syms[2]), syms[3:end]...)
    end
end
function composite_symmetries(sym::magptsym; order::Integer=1)::magptsym
    if order == 1
        return sym
    elseif order > 1
        return composite_symmetries([sym for i = 1:order]...)
    else
        error("Error occurs at composite_symmetries: order should be a positive integer")
    end
end
function composite_symmetries(syms::magptsym...)::magptsym
    if length(syms) == 1
        return syms[1]
    elseif length(syms) == 2
        nsym = (r=syms[1].r*syms[2].r, s=syms[1].s*syms[2].s, T=syms[1].T*syms[2].T)
        return nsym
    else
        return composite_symmetries(composite_symmetries(syms[1], syms[2]), syms[3:end]...)
    end
end


"""
is_symmetry_unitary(sym::spinptsym) -> Bool
is_symmetry_unitary(sym::magptsym) -> Bool
"""
function is_symmetry_unitary(sym::spinptsym)::Bool
    if det(sym.s) > 0
        return true
    else
        return false
    end
end
function is_symmetry_unitary(sym::magptsym)::Bool
    if sym.T == 1
        return true
    else
        return false
    end
end

"""
change_symmetry_from_spin_basis_to_magnon_basis(sym::spg_symtype; SO2::Int64) -> spg_magnonbasis_symtype
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
function change_symmetry_from_spin_basis_to_magnon_basis(sym::spinptsym; SO2::Vector{<:Real}, tor=1e-5)::magptsym
    #=get rotation matrix=#
    r = sqrt(sum(SO2.^2))
    d = sqrt(sum(SO2[1:2].^2))
    theta = acos(SO2[3]/r)
    if d < tor
        phi = 0.0
    elseif SO2[2] > -tor
        phi = acos(SO2[1]/d)
    else
        phi = 2*pi-acos(SO2[1]/d)
    end
    R = Matrix{Float64}([
         cos(phi)*cos(theta) -sin(phi) cos(phi)*sin(theta)
         sin(phi)*cos(theta)  cos(phi) sin(phi)*sin(theta)
        -sin(theta)           0.0      cos(theta)
    ])

    s_global = deepcopy(sym.s)
    s_local = R*s_global*inv(R)
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
    nsym = (r=sym.r, s=ns, T=T)
    return nsym
end

"""
get_multiplication_table_of_group
"""
function get_multiplication_table_of_group(syms::Union{Vector{spinptsym},Vector{magptsym}})::Matrix{Int64}
    if !isa_group(syms) error("Error occurs at get_multiplication_table_of_group: symmetry operations don't form a group") end
    order = length(syms)
    mul_table = zeros(Int64, order, order)
    for i = 1:order
        for j = 1:order
            nsym = composite_symmetries(syms[i], syms[j])
            idx = findall(sym->is_symmetry_same(sym, nsym), syms)[1]
            mul_table[i, j] = idx
        end
    end
    return mul_table
end

"""
isa_group(syms::Vector{T}) -> Bool where {T<:Union{ptsym, spasym, spinptsym, spinspasym}}
check if a set forms a group
"""
function isa_group(syms0::Vector{T})::Bool where {T<:Union{spinptsym, magptsym}}
    syms = deepcopy(syms0)
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
end #module SpinPointSymmetry
