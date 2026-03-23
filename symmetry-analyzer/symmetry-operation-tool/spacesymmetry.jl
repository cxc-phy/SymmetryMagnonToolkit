try
    typeof(MathTool)
catch UndefVarError
    include("../math-tool/mathtool.jl")
end
try
    typeof(PointSymmetry)
catch UndefVarError
    include("./pointsymmetry.jl")
end
######################################################################################################################
####################################################SpaceSymmetry#####################################################
######################################################################################################################
module SpaceSymmetry
using LinearAlgebra
using ..PointSymmetry
using ..MathTool

export spasym
export get_symmetry_in_kspace
export get_xyz_form_of_symmetry
export get_notation_of_symmetry
export is_Identity
export is_Identity_contained
export is_symmetry_contained
export is_symmetry_same
export get_inverse_of_symmetry
export composite_symmetries
export isa_pointsymmetry
export isa_group

spasym = NamedTuple{(:r, :t), Tuple{Matrix{Int64}, Vector{Float64}}}

"""
get_symmetry_in_kspace(rsym::spasym)
"""
function get_symmetry_in_kspace(rsym::spasym)
    ksym = transpose(inv(rsym.r))
    return ksym
end

"""
get_xyz_form_of_symmetry(sym::spasym; tor=1e-3)
"""
function get_xyz_form_of_symmetry(sym::spasym; tor=1e-3)
    ptform = get_xyz_form_of_point_symmetry(sym.r)
    form = split(ptform[1:end-4], ",")
    if length(form) != 3 error("get_notation_of_symmetry") end
    for idx in [1, 2, 3]
        if sym.t[idx] > tor
            form[idx] = string(form[idx], "+", sym.t[idx])
        elseif -sym.t[idx] > tor
            form[idx] = string(form[idx], sym.t[idx])
        end
    end
    form = join(form, ",")
    return form
end

"""
get_notation_of_symmetry(sym::spasym)
"""
function get_notation_of_symmetry(sym::spasym)
    R = get_notation_of_point_symmetry(sym.r)
    t = join(sym.t, " ")
    notation = string("{", R, "|", t, "}")
    return notation
end

"""
is_Identity(sym::spasym) -> Bool
"""
function is_Identity(sym::spasym)::Bool
    Isym = (r=[1 0 0;0 1 0;0 0 1], t=[0.0,0.0,0.0])
    return is_symmetry_same(Isym, sym, modtrans=false)
end

"""
is_Identity_contained(syms::Vector{spasym}) -> Bool
"""
function is_Identity_contained(syms::Vector{spasym})::Bool
    Isym = (r=[1 0 0;0 1 0;0 0 1], t=[0.0,0.0,0.0])
    return is_symmetry_contained(Isym, syms, modtrans=false)
end

"""
is_symmetry_contained(sym::spasym, syms::Vector{spasym}; modtrans::Bool) -> Bool
"""
function is_symmetry_contained(sym::spasym, syms::Vector{spasym}; modtrans::Bool)::Bool
    for ele in syms
        if is_symmetry_same(ele, sym, modtrans=modtrans) return true end
    end
    return false
end

"""
is_symmetry_same(sym1::spasym, sym2::spasym; modtrans::Bool, tor=1e-3) -> Bool
"""
function is_symmetry_same(sym1::spasym, sym2::spasym; modtrans::Bool, tor=1e-3)::Bool
    if modtrans
        relative_t = mod.(round.(sym1.t-sym2.t, digits=2), 1)
    else
        relative_t = round.(sym1.t-sym2.t, digits=2)
    end
    if sym1.r == sym2.r && all(abs.(relative_t) .< tor)
        return true
    else
        return false
    end
end

"""
get_inverse_of_symmetry(sym::spasym; modtrans::Bool) -> spasym
"""
function get_inverse_of_symmetry(sym::spasym; modtrans::Bool)::spasym
    r_inv = get_inverse_of_matrix(sym.r)
    if modtrans
        sym_inv = (r=r_inv, t=mod.(-r_inv*sym.t, 1))
    else
        sym_inv = (r=r_inv, t=-r_inv*sym.t)
    end
    return sym_inv
end

"""
composite_symmetries(syms::spasym; order::Integer, modtrans::Bool) -> spasym
composite_symmetries(syms::spasym...; modtrans::Bool) -> spasym
"""
function composite_symmetries(sym::spasym; order::Integer, modtrans::Bool)::spasym
    if order == 1
        return sym
    elseif order > 1
        return composite_symmetries([sym for i = 1:order]..., modtrans=modtrans)
    else
        error("Error occurs at composite_symmetries: order should be a positive integer")
    end
end
function composite_symmetries(syms::spasym...; modtrans::Bool)::spasym
    if length(syms) == 1
        return syms[1]
    elseif length(syms) == 2
        if modtrans
            return (r=syms[1].r*syms[2].r, t=mod.(syms[1].t+syms[1].r*syms[2].t,1))
        else
            return (r=syms[1].r*syms[2].r, t=syms[1].t+syms[1].r*syms[2].t)
        end
    else
        return composite_symmetries(composite_symmetries(syms[1], syms[2], modtrans = modtrans), syms[3:end]..., modtrans = modtrans)
    end
end


"""
isa_pointsymmetry(sym::spasym) -> Bool
"""
function isa_pointsymmetry(sym::spasym)::Bool
    order = get_symmetry_order_of_point_symmetry(sym.r)
    nsym = composite_symmetries(sym, order=order, modtrans=false)
    if all(abs.(round.(nsym.t, digits=3)) .< 1e-3)
        return true
    else
        return false
    end
end

"""
isa_group(syms::Vector{spasym}) -> Bool
check if a set forms a group
"""
function isa_group(syms::Vector{spasym})::Bool
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

end#(module SpaceSymmetry)
