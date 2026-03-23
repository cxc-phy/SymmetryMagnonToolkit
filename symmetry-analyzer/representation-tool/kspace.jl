try
    typeof(MathTool)
catch UndefVarError
    include("../math-tool/mathtool.jl")
end
######################################################################################################################
########################################################kSpace########################################################
######################################################################################################################
module kSpace
using LinearAlgebra
using ..MathTool

export transform_k
export get_kxkykz_form_of_symmetry
export get_k_invariant_under_group
export get_k_invariant_under_symmetry_operation
export get_intersection_of_k
export is_k_invariant_under_group
export is_k_invariant_under_symmetry_operation
export get_k_geometric_type
export is_k_equivalent_under_group
export is_k_equivalent
export is_k1_contained_in_k2
export is_k_in_klist

ksymtype = Matrix{Float64}
ktype = Matrix{Float64}

"""
transform_k(k::ktype, Tmat::Matrix{<:Number}) -> ktype
    (rx, ry, rz)   ->  (rx, ry, rz) * Tmat
    (x, y, z)^T    ->  Tmat^-1 * (x, y, z)^T
    (kx, ky, kz)^T ->  Tmat^T * (kx, ky, kz)^T
"""
function transform_k(k::ktype, Tmat::Matrix{<:Number})::ktype
    return transpose(Tmat[1:3, 1:3])*k
end

"""
get_kxkykz_form_of_symmetry(ksym::ksymtype; tor=1e-3)
"""
function get_kxkykz_form_of_symmetry(ksym::ksymtype; tor=1e-3)
    form = ["", "", ""]
    kxyz = ["kx", "ky", "kz"]
    ksym = round.(ksym, digits=5)
    for row in 1:3
        for idx in 1:3
            if ksym[row, idx] > 0
                if length(strip(form[row])) != 0 form[row] = string(form[row], "+") end
                if !isinteger(ksym[row, idx]) form[row] = string(form[row], ksym[row, idx]) end
                form[row] = string(form[row], kxyz[idx])
            elseif ksym[row, idx] < 0
                form[row] = string(form[row], "-")
                if !isinteger(ksym[row, idx]) form[row] = string(form[row], -ksym[row, idx]) end
                form[row] = string(form[row], kxyz[idx])
            end
        end
    end
    form = join(form, ", ")
    return form
end

"""
get_k_invariant_under_group(ksyms::Vector{ksymtype}) -> Vector{ktype}
"""
function get_k_invariant_under_group(ksyms::Vector{ksymtype})::Vector{ktype}
    k_list = [Float64.([1 0 0 0;0 1 0 0;0 0 1 0])]
    for ksym in ksyms
        nk_list = []
        for k1 in k_list, k2 in get_k_invariant_under_symmetry_operation(ksym)
            for k in get_intersection_of_k(k1, k2)
                if all([!is_k_equivalent(k, nk, trans=true) for nk in nk_list])
                    push!(nk_list, k)
                end
            end
        end
        k_list = nk_list[:]
    end
    return k_list
end

"""
get_k_invariant_under_symmetry_operation(ksym::ksymtype) -> Vector{ktype}
"""
function get_k_invariant_under_symmetry_operation(ksym::ksymtype)::Vector{ktype}
    k_init = Float64.([1 0 0 0;0 1 0 0;0 0 1 0])
    kconstrs = k_init-ksym*k_init
    kconstrs_list = get_constraints_list_of_symmetry(kconstrs)

    return add_constraints_to_k(k_init, kconstrs_list)
end

"""
get_intersection_of_k(k1, k2; tor=1e-3)
(0,0)(0,1)(0,2)(0,3)
(1,1)(1,2)(1,3)
(2,2)(2,3)
(3,3)
"""
function get_intersection_of_k(k1, k2; tor=1e-3)
    k1dim = get_k_geometric_type(k1)
    k2dim = get_k_geometric_type(k2)
    if k1dim > k2dim
        kdim = k1dim
        k1dim = k2dim
        k2dim = kdim
        k = k1[:,:]
        k1 = k2[:,:]
        k2 = k[:,:]
    end
    if k1dim == k2dim
        if is_k_equivalent(k1, k2, trans=true)
            @goto return_k1
        else
            @goto get_intersection_of_k1_and_k2
        end
    else
        if is_k1_contained_in_k2(k1=k1, k2=k2)
            @goto return_k1
        elseif k1dim == 0
            @goto no_result
        else
            @goto get_intersection_of_k1_and_k2
        end
    end

    @label return_k1
    return [k1]

    @label get_intersection_of_k1_and_k2
    kconstrs = k1 - k2
    kconstrs_list = get_constraints_list_of_symmetry(kconstrs)
    if kconstrs_list == nothing @goto no_result end
    return add_constraints_to_k(k1, kconstrs_list)

    @label no_result
    return []
end

"""
add_constraints_to_k(k, kconstrs_list; tor=1e-3) -> Vector{ktype}
"""
function add_constraints_to_k(k, kconstrs_list; tor=1e-3)::Vector{ktype}
    k_list = Vector{ktype}()
    for kconstrs in kconstrs_list
        nk = k[:, :]
        for var in [3, 2, 1]
            if all(abs.(round.(kconstrs[var, :], digits=3)) .< tor) continue end
            for line in [1, 2, 3]
                nk[line, :] .-= nk[line, var]/kconstrs[var, var]*kconstrs[var, :]
            end
        end
        nk[:, 4] = mod.(nk[:, 4], 1)
        push!(k_list, nk)
    end
    return k_list
end

"""
get_constraints_list_of_symmetry(kconstrs; tor=1e-3)
"""
function get_constraints_list_of_symmetry(kconstrs; tor=1e-3)
    #=
    kconstrs:store constraints of k
    size(kconstrs)=(3,4) ->
        kconstrs[1,1]*u+kconstrs[1,2]*v+kconstrs[1,3]*w+kconstrs[1,4]=0
        kconstrs[2,1]*u+kconstrs[2,2]*v+kconstrs[2,3]*w+kconstrs[1,4]=0
        kconstrs[3,1]*u+kconstrs[3,2]*v+kconstrs[3,3]*w+kconstrs[1,4]=0
    nkconstrs_list:store several reduced constraints of k
    nkconstrs:store one reduced constraints of k
    size(nkconstrs)=(3,4) ->
        kconstrs[1,1]*u+kconstrs[1,4]=0
        kconstrs[2,1]*u+kconstrs[2,2]*v+kconstrs[1,4]=0
        kconstrs[3,1]*u+kconstrs[3,2]*v+kconstrs[3,3]*w+kconstrs[1,4]=0
    used_row: store the rowth row that has already been reduced
    =#
    for row in [1, 2, 3]
        if all(abs.(round.(kconstrs[row, 1:3], digits=3)) .< tor) && abs(round(kconstrs[row, 4], digits=3)) > tor
            @goto no_result
        end
    end
    nkconstrs_list = []
    used_row = []
    for var in [3, 2, 1]
        for row in [3, 2, 1]
            if abs(round(kconstrs[row, var], digits=3)) < tor || row in used_row
                #=
                condition: the rowth row gives no constraints on varth variable or is already used
                statement: continue
                =#
                continue
            end
            if all([i==var ? 0 : abs(round(kconstrs[row, i], digits=3)) for i in 1:3] .< tor)
                #=
                condition: the rowth row in kconstrs can give exact values of the varth variable
                statement:
                #1.put the rowth row into used_row
                #2.get exact values
                #3.go through other rows not in used_row, check if each value is compatible with them
                    if true,
                        (1)store value in kvarsols
                        (2)if the row gives the same result as current one, put it into used_row
                    if false, do nothing
                #4.append rows in used_idx into used_row
                kvarsols: store solution for the varth variable
                used_idx: store the idxth row giving the same constraints as the rowth row
                =#
                push!(used_row, row)
                kvarsols = []
                used_idx = []
                for sol in  get_exact_solution_of_k_component(var, kconstrs[row, :])
                    kvarsol = [[i == var ? -1 : 0 for i = 1:3];sol]
                    for idx in [3, 2, 1]
                        if abs(round(kconstrs[idx, var], digits=3)) < tor || idx in used_row continue end
                        constr = kconstrs[idx, :] - kconstrs[idx, var]/kvarsol[var]*kvarsol[:]
                        if any(abs.(round.(constr[1:3], digits=3)) .> tor)
                            #kconstrs[idx, :] .-= kconstrs[idx, var]/kvarsol[var]*kvarsol[:]
                            continue
                        end
                        if !(idx in used_idx) push!(used_idx, idx) end
                        if abs(mod(round(constr[4], digits=3), 1)) < tor continue end
                        @goto next_sol
                    end
                    push!(kvarsols, kvarsol)
                    @label next_sol
                end
                if kvarsols == [] @goto no_result end
                append!(used_row, used_idx)
            else
                #=
                condition: the rowth row in kconstrs gives constraint on the varth variable
                statement:
                #1.check if the coefficent of the varth variable (kconstrs[row, var]) is +1/-1
                    if true, store constraint in kvarsols
                    if false and exist another row having true result, do nothing
                    if false and exist no other row having true result, store constraint in kvarsols
                #2.if constraint is stored, eliminate the varth variable in the remaining constraints
                *the purpose is to avoid division due to the constraints is actually mod 1
                kvarsols: store solution for the varth variable
                =#
                if abs(round(abs(kconstrs[row, var])-1, digits=3)) > tor
                    @goto next_row
                    for idx = 1:row-1
                        if abs(round(abs(kconstrs[idx, var])-1, digits=3)) < tor @goto next_row end
                    end
                end
                push!(used_row, row)
                kvarsols = []
                push!(kvarsols, kconstrs[row, :]*-kconstrs[row, var]/abs(kconstrs[row, var]))
                for idx in [3, 2, 1]
                    if abs(round(kconstrs[idx, var], digits=3)) < tor || idx in used_row continue end
                    kconstrs[idx, :] .-= kconstrs[idx, var]/kconstrs[row, var]*kconstrs[row, :]
                end
            end
            #=

            add kvarsols to nkconstrs_list
            =#
            if nkconstrs_list == []
                nkconstrs_list = [kvarsol for kvarsol in kvarsols]
            else
                nkconstrs_list = [[kvarsol; nkconstrs] for kvarsol in kvarsols for nkconstrs in nkconstrs_list]
            end
            break
            @label next_row
        end
        #=the varth variable has no exact values or constraints=#
        if nkconstrs_list == []
            nkconstrs_list = [Float64.([0, 0, 0, 0])]
        elseif length(nkconstrs_list[1]) == 4*(3-var)
            nkconstrs_list = [[Float64.([0, 0, 0, 0]); nkconstrs] for nkconstrs in nkconstrs_list]
        elseif length(nkconstrs_list[1]) != 4*(4-var)
            error("Error occurs at get_constraints_list_of_symmetry")
        end
    end
    nkconstrs_list = [[ele[4*(i-1)+j] for i=1:3, j=1:4] for ele in nkconstrs_list]
    for var in [1, 2]
        for i in eachindex(nkconstrs_list)
            nkconstrs = nkconstrs_list[i]
            if abs(round(nkconstrs[var, var], digits=3)) < tor continue end
            nkconstrs = [ idx in var+1:3 ? nkconstrs[idx, :] .-= nkconstrs[idx, var]/nkconstrs[var, var]*nkconstrs[var, :] : nkconstrs[idx, :] for idx = 1:3]
            nkconstrs_list[i] = [nkconstrs[i][j] for i=1:3, j=1:4]
        end
    end
    return nkconstrs_list

    @label no_result
    return nothing
end

"""
get_exact_solution_of_k_component(var, kconstr; tor=1e-3)
"""
function get_exact_solution_of_k_component(var, kconstr; tor=1e-3)
    exactsol = []
    for n = -6:6
        ans = -(n+kconstr[4])/kconstr[var]
        if round(ans-0.5, digits=3) > tor || round(-0.5-ans, digits=3) > -tor continue end
        for sol in exactsol
            if abs(round(ans-sol, digits=3)) < tor @goto next_integer end
        end
        push!(exactsol, ans)
        @label next_integer
    end
    return exactsol
end

"""
is_k_invariant_under_group(k::ktype, ksyms::Vector{ksymtype}) -> Bool
    G = true: k and nk can differ by a reciprocal lattice vector
"""
function is_k_invariant_under_group(k::ktype, ksyms::Vector{ksymtype}; G::Bool=true)::Bool
    for ksym in ksyms
        if !is_k_invariant_under_symmetry_operation(k, ksym, G=G) return false end
    end
    return true
end

"""
is_k_invariant_under_symmetry_operation(k::ktype, ksym::ksymtype) -> Bool
#check if k is invariant under given symmetry
    G = true: k and nk can differ by a reciprocal lattice vector
#methods
    check if new k is exactly the same as k (each kpoint should transform to itself)
"""
function is_k_invariant_under_symmetry_operation(k::ktype, ksym::ksymtype; G::Bool=true)::Bool
    nk = ksym * k
    return is_k_equivalent(k, nk, trans=false, G=G)
end

"""
get_k_geometric_type(k::ktype; tor=1e-3) -> Integer
#get geometric type of k -> 0: kpoint; 1: kline; 2: kplane; 3: kspace
#method: find number of free variables
#conditions for two dependent variables k[i,1:3] and k[j,1:3]
    1. k[i,m]/k[j,m]!=NaN
    2. k[i,1]/k[j,1]=k[i,2]/k[j,2]=k[i,3]/k[j,3]
"""
function get_k_geometric_type(k::ktype; tor=1e-3)::Integer
    if all(abs.(round.(k[:, 1:3], digits=3)) .< tor) return 0 end
    var_idx = []
    for i = 1:3
        if all(abs.(round.(k[i, 1:3], digits=3)) .< tor) continue end
        if var_idx == []
            push!(var_idx, i)
        else
            independent = true
            for idx in var_idx
                kratio = [k[i, j]/k[idx, j] for j=1:3 if !isequal(k[i, j]/k[idx, j], NaN)]
                if length(kratio) < 2
                    independent = false
                    break
                elseif Inf in kratio
                    continue
                else
                    for i = 1:length(kratio)-1
                        if all(abs.(round.(kratio[i]-kratio[i+1], digits=3)) .< tor)
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
    return length(var_idx)
end

"""
is_k_equivalent_under_group(k1::ktype, k2::ktype, ksyms::Vector{ksymtype}) -> Bool
"""
function is_k_equivalent_under_group(k1::ktype, k2::ktype, ksyms::Vector{ksymtype}; trans::Bool=true, G::Bool=true)::Bool
    k1type = get_k_geometric_type(k1)
    k2type = get_k_geometric_type(k2)
    if k1type != k2type return false end
    for ksym in ksyms
        nk1 = ksym * k1
        if is_k_equivalent(nk1, k2, trans=trans, G=G) return true end
    end
    return false
end


"""
is_k_equivalent(k1::ktype, k2::ktype; trans::Bool) -> Bool
#check if two k are equivalent
#methods
    1.check if Reals of free variable are same
    2.go to corresponding function is_(kpoint/kline/kplane/kspace)_equivalent
"""
function is_k_equivalent(k1::ktype, k2::ktype; trans::Bool, G::Bool=true)::Bool
    k1type = get_k_geometric_type(k1)
    k2type = get_k_geometric_type(k2)
    if k1type != k2type return false end
    if k1type == 0 return is_kpoint_equivalent(k1, k2, G=G) end
    if k1type == 1 return is_kline_equivalent(k1, k2, trans=trans, G=G) end
    if k1type == 2 return is_kplane_equivalent(k1, k2, trans=trans, G=G) end
    if k1type == 3 return is_kspace_equivalent(k1, k2, trans=trans) end
end

"""
is_kpoint_equivalent(kp1::ktype, kp2::ktype; tor=1e-3) -> Bool
is_kpoint_equivalent(kp1::Vector{Float64}, kp2::Vector{Float64}; tor=1e-3) -> Bool
#check if two kpoints are equivelent
    1. G = true, kp1 and kp2 can differ by a reciprocal lattice vector
#size(kp) = (3,) -> (
    kx = kp[1],
    ky = kp[2],
    kz = kp[3]
)
"""
function is_kpoint_equivalent(kpoint1::ktype, kpoint2::ktype; G::Bool=true, tor=1e-3)::Bool
    kp1 = kpoint1[:, 4]
    kp2 = kpoint2[:, 4]
    return is_kpoint_equivalent(kp1, kp2, G=G)
end
function is_kpoint_equivalent(kpoint1::Vector{Float64}, kpoint2::Vector{Float64}; G=true, tor=1e-3)::Bool
    kp1 = deepcopy(kpoint1)
    kp2 = deepcopy(kpoint2)
    if G
        kp1 = mod.(kp1, 1)
        kp2 = mod.(kp2, 1)
    end
    if all(abs.(round.(kp1-kp2,digits=3)).<[tor tor tor])
        return true
    else
        return false
    end
end

"""
is_kline_equivalent(kl1::ktype, kl2::ktype; trans, tor=1e-3) -> Bool
#check if two klines are equivelent
    1. trans = true: given certain (u,v,w), each point on kl1 can be different from kl2 up to an arbitrary global translation along the kline
    2. trans = false: given certain (u,v,w), each point on kl1 should be exactly the same as kl2
    3. G = true, kl1 and kl2 can differ by a reciprocal lattice vector
#conditions
    1. have same direction vectors
    2.(1) trans = true: given (u,v,w), have one same kpoint up to an arbitrary global translation along the kline
    2.(2) trans = false: given (u,v,w), have one same kpoint
#size(kl) = (3, 4) -> (
    kx = kl[1,1]*u+kl[1,2]*v+kl[1,3]*w+kl[1,4],
    ky = kl[2,1]*u+kl[2,2]*v+kl[2,3]*w+kl[2,4],
    kz = kl[3,1]*u+kl[3,2]*v+kl[3,3]*w+kl[3,4]
)
"""
function is_kline_equivalent(kline1::ktype, kline2::ktype; trans::Bool, G::Bool=true, tor=1e-3)::Bool
    kl1 = deepcopy(kline1)
    kl2 = deepcopy(kline2)
    if G
        kl1[:, 4] = mod.(kl1[:, 4], 1)
        kl2[:, 4] = mod.(kl2[:, 4], 1)
    end

    #=check if direction vectors are equivalent=#
    dirvec1 = MathTool.get_direction_vector_of_line(kl1)
    dirvec2 = MathTool.get_direction_vector_of_line(kl2)
    if dirvec1 == nothing || dirvec2 == nothing error("Error occurs at is_kline_equivalent: can not find direction vector") end
    if any(abs.(round.(dirvec1-dirvec2,digits=3)).>[tor tor tor]) && any(abs.(round.(dirvec1+dirvec2,digits=3)).>[tor tor tor]) return false end

    #=
    check if having one same kpoint
    1.get kpoints at (u=0, v=0, w=0)
    2.(1)if trans = false, check if two kpoints are equivalent
    2.(2)if trans = true, check if direction vector between two kpoints is paralleled to the direction vector
    =#
    uvw = [0.1, 0.2, 0.3, 1]
    kp1 = kl1*uvw
    kp2 = kl2*uvw
    if trans
        dirvec12 = kp1 - kp2
        if any(abs.(round.(cross(dirvec12, dirvec1), digits=3)) .>[tor tor tor])
            return false
        else
            return true
        end
    else
        return is_kpoint_equivalent(kp1, kp2, G=G)
    end

end

"""
is_kplane_equivalent(kp1::ktype, kp2::ktype; trans, tor=1e-3) -> Bool
#check if two kplanes are equivelent
    1. trans = true: given certain (u,v,w), each point on kp1 can be different from kp2 up to an arbitrary global translation along the kplane
    2. trans = false: given certain (u,v,w), each point on kp1 should be exactly the same as kp2
    3. G = true, kpl1 and kpl2 can differ by a reciprocal lattice vector
#conditions
    1. have same normal vectors
    2.(1) trans = true: given (u,v,w), have one same kpoint up to an arbitrary global translation along the kplane
    2.(2) trans = false: given (u,v,w), have one same kpoint
#size(kp) = (3, 4) -> (
    kx = kp[1,1]*u+kp[1,2]*v+kp[1,3]*w+kp[1,4],
    ky = kp[2,1]*u+kp[2,2]*v+kp[2,3]*w+kp[2,4],
    kz = kp[3,1]*u+kp[3,2]*v+kp[3,3]*w+kp[3,4]
)
"""
function is_kplane_equivalent(kplane1::ktype, kplane2::ktype; trans::Bool, G::Bool=true, tor=1e-3)::Bool
    kpl1 = deepcopy(kplane1)
    kpl2 = deepcopy(kplane2)
    if G
        kpl1[:, 4] = mod.(kpl1[:, 4], 1)
        kpl2[:, 4] = mod.(kpl2[:, 4], 1)
    end

    #=check if normal vectors are equivalent=#
    normvec1 = MathTool.get_normal_vector_of_plane(kpl1)
    normvec2 = MathTool.get_normal_vector_of_plane(kpl2)
    if normvec1 == nothing || normvec2 == nothing error("Error occurs at is_kplane_equivalent: can not find normal vector") end
    if any(abs.(round.(normvec1-normvec2,digits=3)).>[tor tor tor]) && any(abs.(round.(normvec1+normvec2,digits=3)).>[tor tor tor]) return false end

    #=
    check if having one same kpoint
    1.get kpoints at (u=0, v=0, w=0)
    2.(1)if trans = false, check if two kpoints are equivalent
    2.(2)if trans = true, check if direction vector between two kpoints is normal to the normal vector
    =#
    uvw = [0.1, 0.2, 0.3, 1]
    kp1 = kpl1*uvw
    kp2 = kpl2*uvw
    #kp1start = mod.(kp1[:, 4], 1)
    #kp2start = mod.(kp2[:, 4], 1)
    if trans
        dirvec12 = kp1 - kp2
        if abs(round(dot(dirvec12, normvec1), digits=3)) > tor
            return false
        else
            return true
        end
    else
        return is_kpoint_equivalent(kp1, kp2, G=G)
    end

end


"""
is_kspace_equivalent(ksp1::ktype, ksp2::ktype; trans, tor=1e-3) -> Bool
#check if two kspace are equivelent
    1. trans = true: given certain (u,v,w), each point on kl1 can be different from kl2 up to an arbitrary global translation(naturally satisfied)
    2. trans = false: given certain (u,v,w), each point on kl1 should be exactly the same as kl2
#conditions
    1.(1) trans = true: given (u,v,w), have one same kpoint up to an arbitrary global translation(naturally satisfied)
    1.(2) trans = false: given (u,v,w), have one same kpoint
#size(kl) = (3, 4) -> (
    kx = kl[1,1]*u+kl[1,2]*v+kl[1,3]*w+kl[1,4],
    ky = kl[2,1]*u+kl[2,2]*v+kl[2,3]*w+kl[2,4],
    kz = kl[3,1]*u+kl[3,2]*v+kl[3,3]*w+kl[3,4]
)
"""
function is_kspace_equivalent(kspace1::ktype, kspace2::ktype; trans::Bool, tor=1e-3)::Bool
    ksp1 = deepcopy(kspace1)
    ksp2 = deepcopy(kspace2)
    ksp1[:, 4] = mod.(ksp1[:, 4], 1)
    ksp2[:, 4] = mod.(ksp2[:, 4], 1)

    if trans
        return true
    else
        uvwlist = [
            [0.0, 0.0, 0.0, 1],
            [0.0, 0.0, 0.1, 1],
            [0.0, 0.1, 0.0, 1],
            [0.0, 0.1, 0.1, 1],
            [0.1, 0.0, 0.0, 1],
            [0.1, 0.0, 0.1, 1],
            [0.1, 0.1, 0.0, 1],
            [0.1, 0.1, 0.1, 1]
        ]
        for uvw in uvwlist
            kp1 = ksp1*uvw
            kp2 = ksp2*uvw
            if !is_kpoint_equivalent(kp1, kp2) return false end
        end
        return true
    end

end

"""
is_k1_contained_in_k2(;k1, k2) -> Bool
(0,1)(0,2)(0,3)
(1,2)(1,3)
(2,3)
"""
function is_k1_contained_in_k2(;k1, k2)::Bool
    k1dim = get_k_geometric_type(k1)
    k2dim = get_k_geometric_type(k2)
    if k1dim >= k2dim
        error("dimension of geometric type of k1 should be smaller than k2")
    end
    if k2dim == 3 return true end
    if k1dim == 0
        return is_kpoint_contained_in_k(k1, k2)
    elseif k1dim == 1
        return is_kline_contained_in_k(k1, k2)
    end
end

"""
is_kpoint_contained_in_k(kp, k; tor=1e-3) -> Bool
(0,1)(0,2)
"""
function is_kpoint_contained_in_k(kp, k; tor=1e-3)::Bool
    for line in [1,2,3]
        if all(abs.(round.(k[line, 1:3], digits=3)) .< tor)
            if abs(round(mod(k[line, 4]-kp[line, 4], 1), digits=3)) > tor
                return false
            end
        elseif abs(round(sum(k[line, 1:3])-k[line, line], digits=3)) > tor
            if abs(round(mod(sum(kp[i, 4]*k[line, i] for i = 1:3)+k[line, 4]-kp[line, 4], 1), digits=3)) > tor
                return false
            end
        end
    end
    return true
end

"""
is_kline_contained_in_k(kl, k; tor=1e-3) -> Bool
(1,2)
"""
function is_kline_contained_in_k(kl, k; tor=1e-3)::Bool
    dirvec = get_direction_vector_of_line(kl)
    normvec = get_normal_vector_of_plane(k)
    if abs(round(sum(dirvec.*normvec), digits=3)) > tor
        return false
    end
    kp = [zeros(Float64, 3, 3);;kl[:, 4]]
    return is_kpoint_contained_in_k(kp, k)
end

"""
is_k_in_klist(k, klist)::Bool
"""
function is_k_in_klist(k, klist)::Bool
    for q in klist
        if is_k_equivalent(k, q, trans=true) return true end
    end
    return false
end

end #(module kSpace)
