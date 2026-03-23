try
    typeof(MathTool)
catch UndefVarError
    include("../math-tool/mathtool.jl")
end
try
    typeof(Symmetry)
catch UndefVarError
    include("../symmetry-operation-tool/symmetry.jl")
end
######################################################################################################################
######################################################RealSpace#######################################################
######################################################################################################################
module RealSpace
using LinearAlgebra
using ..MathTool
import ..Symmetry as sym

export rtype
export get_fixedpoints_of_symmetry_operation
export get_r_invariant_under_symmetry_operation
export get_intersection_of_r
export is_r1_contained_in_r2
export is_r_same
export get_r_geometric_type

rtype = Matrix{Float64}
rsymtype = NamedTuple{(:r, :t), Tuple{Matrix{Int64}, Vector{Float64}}}
#=
get_intersection_of_r(r1::rtype, r2::rtype; tor=1e-3) -> rtype
is_r1_contained_in_r2(;r1, r2) -> Bool
    is_rpoint_contained_in_r(rp, r; tor=1e-3) -> Bool
    is_rline_contained_in_r(rl, r; tor=1e-3) -> Bool
is_r_same(r1::rtype, r2::rtype; tor=1e-3) -> Bool
get_r_geometric_type(r::rtype; tor=1e-3) -> Int64
=#

function get_fixedpoints_of_symmetry_operation(rsym::rsymtype; tor=1e-3)
    t_list = [-2.0, -1.0, 0.0, 1.0, 2.0]
    fixedpoints = []
    for tx in t_list
        for ty in t_list
            for tz in t_list
                nrsym = (r=rsym.r, t=rsym.t+[tx,ty,tz])
                if !sym.isa_pointsymmetry(nrsym) continue end
                fixedpoint = get_r_invariant_under_symmetry_operation(nrsym)
                for (sym, fixed) in fixedpoints
                    if all(abs.(round.(fixed-fixedpoint, digits=3)) .< tor) @goto next end
                end
                push!(fixedpoints, (nrsym, fixedpoint))
                @label next
            end
        end
    end
    return fixedpoints
end

"""
get_r_invariant_under_symmetry_operation(rsym::rsymtype) -> Matrix{Float64}
"""
function get_r_invariant_under_symmetry_operation(rsym::rsymtype; tor=1e-3)
    A = rsym.r .- Matrix{Int64}(I, 3, 3)
    b = -rsym.t
    fixedpoint = nothing
    constrs = nothing
    try
        x = A\b #round.(A\b, digits=3)
        fixedpoint = zeros(Float64, 3, 4)
        for i in eachindex(x)
            if isequal(x[i], NaN)
                fixedpoint[i, i] = 1
            elseif isequal(x[i], Inf) || isequal(x[i], -Inf)
                error("not a symmorphic symmetry")
            else
                fixedpoint[i, 4] = x[i]
            end
        end
    catch e
        constrs = trape_mat([A;;-b])
        #constrs = -optimize_constraints([A;;-b])[1]
        fixedpoint = solve_constraint(constrs)
    end

    #check
    if !all(abs.(round.(rsym.r*fixedpoint+[zeros(3,3);;rsym.t]-fixedpoint, digits=3)) .< tor)
        println(000,constrs)
        println(111,rsym)
        println(222,fixedpoint)
        println(333,rsym.r*fixedpoint+[zeros(3,3);;rsym.t]-fixedpoint)
        error()
    end

    return fixedpoint
end

function swap_row(a,i,j)
    m,n=size(a)
    if i>m || j >m
        println(a)
        println(i,j)
        error()
    else
        row = a[i,:]
        a[i,:]=a[j,:]
        a[j,:]=row
    end
    return a
end

function trape_mat(sigma)
    constr = Float64.([0 0 0 0;0 0 0 0;0 0 0 0])
    m,n = size(sigma)
    main_factor = []
    main_col = n-1
    while (main_col > 0 && length(main_factor) < m)
        last_row = m-length(main_factor)
        not_zeros = [i for i in 1:last_row if sigma[i, main_col] != 0]
        if length(not_zeros) != 0
            push!(main_factor, main_col)
            index = not_zeros[end]
            constr[main_col, :] = sigma[index,:]
            if index != last_row sigma = swap_row(sigma, last_row, index) end
            if last_row != 1
                for k in range(last_row-1, 1,step=-1)
                    times = sigma[k, main_col]/sigma[last_row, main_col]
                    sigma[k,:]=sigma[k,:]-times*sigma[last_row,:]
                end
            end
        end
        main_col -=1
    end
    return constr
end

"""
add_constraints_to_k(k, constrs_list; tor=1e-3) -> Vector{ktype}
"""
function solve_constraint(constrs; tor=1e-3)::rtype
    r = Float64.([1 0 0 0;0 1 0 0;0 0 1 0])
    for var in [3, 2, 1]
        if all(abs.(round.(constrs[var, :], digits=3)) .< tor) continue end
        for line in [1, 2, 3]
            r[line, :] .-= r[line, var]/constrs[var, var]*constrs[var, :]
        end
    end
    #r[:, 4] = mod.(r[:, 4], 1)
    return r
end

"""
get_intersection_of_r(r1::rtype, r2::rtype; tor=1e-3) -> rtype
"""
function get_intersection_of_r(r1::rtype, r2::rtype; tor=1e-3)::Union{rtype, Nothing}
    try
        result = round.((r1[:, 1:3]-r2[:, 1:3])\(r2[:, 4]-r1[:, 4]), digits=3)
        intersec = zeros(Float64, 3, 4)
        for i in eachindex(result)
            if isequal(result[i], NaN)
                intersec[i, :] = r1[i, :]
            elseif isequal(result[i], Inf) || isequal(result[i], -Inf)
                return nothing
            else
                intersec[i, 4] = result[i]
            end
        end
        return intersec
    catch e
        if all(abs.(round.((r1[:, 1:3]-r2[:, 1:3]), digits=3)) .< tor) && all(abs.(round.((r2[:, 4]-r1[:, 4]), digits=3)) .< tor) return r1 end
        #println("Warning:get_intersection_of_r")
        #println("    ",r1[:, 1:3]-r2[:, 1:3])
        #println("    ",r2[:, 4]-r1[:, 4])
        return nothing
    end
    #=
    intersec = zeros(Float64, 3, 4)
    for i in eachindex(result)
        if isequal(result[i], NaN)
            intersec[i, :] = r1[i, :]
        elseif isequal(result[i], Inf) || isequal(result[i], -Inf)
            return nothing
        else
            intersec[i, 4] = result[i]
        end
    end
    return intersec
    =#
end

"""
is_r1_contained_in_r2(;r1, r2) -> Bool
(0,1)(0,2)(0,3)
(1,2)(1,3)
(2,3)
"""
function is_r1_contained_in_r2(;r1::rtype, r2::rtype)::Bool
    r1dim = get_r_geometric_type(r1)
    r2dim = get_r_geometric_type(r2)
    if r1dim >= r2dim
        error("dimension of geometric type of r1 should be smaller than r2")
    end
    if r2dim == 3 return true end
    if r1dim == 0
        return is_rpoint_contained_in_r(r1, r2)
    elseif r1dim == 1
        return is_rline_contained_in_r(r1, r2)
    end
end

"""
is_rpoint_contained_in_r(rp, r; tor=1e-3) -> Bool
(0,1)(0,2)
"""
function is_rpoint_contained_in_r(rp, r; tor=1e-3)::Bool
    for line in [1,2,3]
        if all(abs.(round.(r[line, 1:3], digits=3)) .< tor)
            if abs(round(r[line, 4]-rp[line, 4], digits=3)) > tor
                return false
            end
        elseif abs(round(sum(r[line, 1:3])-r[line, line], digits=3)) > tor
            if abs(round(sum(rp[i, 4]*r[line, i] for i = 1:3)+r[line, 4]-rp[line, 4], digits=3)) > tor
                return false
            end
        end
    end
    return true
end

"""
is_rline_contained_in_r(rl, r; tor=1e-3) -> Bool
(1,2)
"""
function is_rline_contained_in_r(rl, r; tor=1e-3)::Bool
    dirvec = get_direction_vector_of_line(rl)
    normvec = get_normal_vector_of_plane(r)
    if abs(round(sum(dirvec.*normvec), digits=3)) > tor
        return false
    end
    rp = [zeros(Float64, 3, 3);;rl[:, 4]]
    return is_rpoint_contained_in_r(rp, r)
end

"""
is_r_same(r1::rtype, r2::rtype; tor=1e-3) -> Bool
"""
function is_r_same(r1::rtype, r2::rtype; modtrans::Bool, tor=1e-3)::Bool
    relative_r = abs.(round.(r1-r2, digits=3))
    if modtrans
        relative_r[:, 4] = mod.(relative_r[:, 4], 1)
    end
    if all(relative_r .< tor)
        return true
    else
        return false
    end
end

"""
get_r_geometric_type(r::rtype; tor=1e-3) -> Int64
"""
function get_r_geometric_type(r::rtype; tor=1e-3)::Int64
    return sum([1 for i in 1:3 if !all(abs.(round.(r[:,i], digits=3)) .< tor)])
    #=if all(abs.(round.(r[:, 1:3], digits=3)) .< tor) return 0 end
    var_idx = []
    for i = 1:3
        if all(abs.(round.(r[i, 1:3], digits=3)) .< tor) continue end
        if var_idx == []
            push!(var_idx, i)
        else
            independent = true
            for idx in var_idx
                rratio = [r[i, j]/r[idx, j] for j=1:3 if !isequal(r[i, j]/r[idx, j], NaN)]
                if length(rratio) < 2
                    independent = false
                    break
                elseif Inf in rratio
                    continue
                else
                    for i = 1:length(rratio)-1
                        if all(abs.(round.(rratio[i]-rratio[i+1], digits=3)) .< tor)
                            independent = false
                            break
                        end
                    end
                end
            end
            if independent push!(var_idx, i) end
        end
    end
    if length(var_idx) == 0 error("Error occurs at get_k_type: k should have at least one variable") end
    return length(var_idx)=#
end

end#(module RealSpace)
