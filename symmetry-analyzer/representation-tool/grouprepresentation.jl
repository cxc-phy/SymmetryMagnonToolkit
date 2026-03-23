try
    typeof(MathTool)
catch UndefVarError
    include("../math-tool/mathtool.jl")
end

module GroupRepresentation
using LinearAlgebra
using ..MathTool

export get_multiplicity_of_irrep
export get_index_of_irrep
export is_irrep_equivalent
export get_min_dimension_of_enhanced_character_table
export get_max_dimension_of_character_table
export get_induced_irreps
export get_similar_transformation_matrix_between_two_reps
export check_transformation_matrix
export get_irrep_basis
export check_irrep_basis
export get_projection_operators
export check_irreps
export get_character_table_of_irreps


"""
get_multiplicity_of_irrep
"""
function get_multiplicity_of_irrep(rep::Vector{Matrix{ComplexF64}}, irrep::Matrix{Matrix{ComplexF64}})::Vector{Int64}
    muls = Int64[]
    for i in axes(irrep, 1)
        mul = get_multiplicity_of_irrep(rep, irrep[i,:])
        push!(muls, mul)
    end
    return muls
end
function get_multiplicity_of_irrep(rep::Vector{Matrix{ComplexF64}}, irrep::Vector{Matrix{ComplexF64}})::Int64
    tr_rep = tr.(rep)
    tr_irrep = tr.(irrep)
    return get_multiplicity_of_irrep(tr_rep, tr_irrep)
end
function get_multiplicity_of_irrep(rep::Vector{Matrix{ComplexF64}}, tr_irrep::Vector{<:Number})::Int64
    tr_rep = tr.(rep)
    return get_multiplicity_of_irrep(tr_rep, tr_irrep)
end
function get_multiplicity_of_irrep(tr_rep::Vector{<:Number}, tr_irrep::Vector{<:Number})::Int64
    mul = complex(0.0, 0.0)
    order = length(tr_rep)
    for i in 1:order
        mul += tr_rep[i]*conj(tr_irrep[i])/order
    end
    #println(mul)
    if !isinteger(round(mul, digits=3)) error("Error occurs at get_multiplicity_of_irrep: multiplicity is not an integer") end
    mul = Int64(round(mul, digits=3))
    return mul
end

"""
get_index_of_irrep(irrep, irreps) -> Int64
"""
function get_index_of_irrep(irrep, irreps)::Int64
    for idx in axes(irreps, 1)
        if is_irrep_equivalent(irrep, irreps[idx, :]) return idx end
    end
    return 0
end
"""
is_irrep_equivalent(irrep1, irrep2) -> Bool
"""
function is_irrep_equivalent(irrep1, irrep2; tor=1e-5)::Bool
    if all(round.(norm.(tr.(irrep1) .- tr.(irrep2)), digits=5) .< tor)
        return true
    else
        return false
    end
end

"""
get_min_dimension_of_enhanced_character_table(table, idx) -> Int64
"""
function get_min_dimension_of_enhanced_character_table(table, idx=0)::Int64
    maxdim = 10
    if typeof(table) <: Matrix{<:Any}
        for i in axes(table, 1)
            if maxdim > Int64(table[i, 1]) maxdim = Int64(table[i, 1]) end
        end
    else
        for i in eachindex(table)
            if idx != 0 && tr(Int64(real(table[i].character[idx]))) > 0 continue end
            if table[i].case == 1 && maxdim > Int64(table[i].character[1])
                maxdim = Int64(table[i].character[1])
            elseif table[i].case != 1 && maxdim > 2*Int64(table[i].character[1])
                maxdim = Int64(2*table[i].character[1])
            end
        end
    end
    return maxdim
end

"""
get_max_dimension_of_character_table(table) -> Int64
"""
function get_max_dimension_of_character_table(table)::Int64
    maxdim = 0
    if typeof(table) <: Matrix{<:Any}
        for i in axes(table, 1)
            if maxdim < Int64(table[i, 1]) maxdim = Int64(table[i, 1]) end
        end
    else
        for i in eachindex(table)
            if table[i].case == 1 && maxdim < Int64(table[i].character[1])
                maxdim = Int64(table[i].character[1])
            elseif table[i].case != 1 && maxdim < 2*Int64(table[i].character[1])
                maxdim = Int64(2*table[i].character[1])
            end
        end
    end
    return maxdim
end

"""
T^-1 * irreps1 * T = irreps2
"""
function get_similar_transformation_matrix_between_two_reps(irreps1, irreps2; tor=1e-3)
    gorder = size(irreps1, 2)
    dim = size(irreps1[1], 1)
    if dim == 1
        Tmatrix = [complex(1.0, 0.0);;]
    else
        if all([all(norm.(round.(irreps1[k]-irreps2[k], digits=5)) .< tor) for k in 1:gorder])
            Tmatrix = complex(1.0, 0.0)*Matrix{Complex{Float64}}(I, dim, dim)
        else
            projectors = get_projection_operators(irreps2, irreps1)
            basis_list = get_irrep_basis(projectors)
            if !check_irrep_basis(basis_list, irreps1, irreps2) error("irrepbasis") end
            Tmatrix = MathTool.get_transformation_matrix(basis_list)
            if !check_transformation_matrix(Tmatrix, irreps1, irreps2) error("Tmatrix") end
        end
    end
    return Tmatrix
end

"""
check_transformation_matrix(Tmatrix, irreps1, irreps2)::Bool
"""
function check_transformation_matrix(Tmatrix, irreps1, irreps2; tor=1e-3)::Bool
    for k in eachindex(irreps1)
        result=composite_matrice(inv(Tmatrix), irreps1[k], Tmatrix)-irreps2[k]
        if !all(norm.(round.(result, digits=5)) .< tor) return false end
    end
    return true
end

"""
get_irrep_basis(projectors; init_basis=nothing, tor=1e-3)
"""
function get_irrep_basis(projectors; init_basis=nothing, tor=1e-3)
    dim = size(projectors[1, 1], 1)
    if init_basis == nothing init_basis = [complex(1.0, 0.0) for i in 1:dim] end
    basis_list = []
    for i in 1:dim
        basis = projectors[i, 1]*init_basis
        if all(norm.(round.(basis, digits=5)) .< tor) error("please change initial basis") end
        push!(basis_list, basis)
    end
    for i in 1:dim
        basisi = projectors[i,i]*basis_list[i]
        if !all(norm.(round.(basisi-basis_list[i], digits=5)) .< tor) error("Error occurs at get_irrep_basis: basis not correct") end
    end
    return basis_list
end

"""
check_irrep_basis(basis_list, reps, irreps; tor=1e-3)::Bool
"""
function check_irrep_basis(basis_list, reps, irreps; tor=1e-3)::Bool
    irrep_dim = size(irreps[1], 1)
    rep_dim = size(reps[1], 1)
    for i in eachindex(irreps)
        rep = reps[i]
        irrep = irreps[i]
        nbasis_list = [rep*basis for basis in basis_list]
        for j in 1:irrep_dim, k in 1:irrep_dim
            if !all(norm.(round.(irrep[j, k] - sum(conj(basis_list[j]) .* nbasis_list[k]), digits=5)).<tor) return false end
        end
    end
    return true
end

"""
get_projection_operators
"""
function get_projection_operators(irreps, reps)
    dim = size(reps[1], 1)
    gorder = length(reps)
    projectors = Matrix{Any}([nothing for i in 1:dim, j in 1:dim])
    for i in 1:dim
        for j in 1:dim
            projector = zeros(Complex{Float64}, dim, dim)
            for idx in eachindex(reps)
                projector .+= conj(irreps[idx][i, j])*reps[idx]
            end
            projectors[i, j] = projector
        end
    end
    projectors = projectors*dim/gorder
    return projectors
end

"""
check_irreps(irreps)::Bool
"""
function check_irreps(irreps; tor=1e-3)::Bool
    #=check number and dimension of irrep=#
    gorder = size(irreps, 2)
    if gorder != sum([size(irreps[i,1], 1)^2 for i in axes(irreps, 1)])
        return false
    end
    #=check wonderful orthogonality theorem for irreps=#
    irrep_num = size(irreps, 1)
    for i in 1:irrep_num
        for j in 1:irrep_num
            dimi = size(irreps[i, 1], 1)
            dimj = size(irreps[j, 1], 1)
            for m in 1:dimi, n in 1:dimi, k in 1:dimj, l in 1:dimj
                result = sum([irreps[i,idx][m,n]* conj(irreps[j,idx][k,l]) for idx in 1:gorder])
                if i == j && m==k && n == l
                    if !all(norm.(round.(result*dimi/gorder-complex(1.0,0.0), digits=5)) .< tor) return false end
                else
                    if !all(norm.(round.(result*dimi/gorder, digits=5)) .< tor) return false end
                end
            end
        end
    end
    character_table, classes = get_character_table_of_irreps(irreps, class=true)
    multi = [length(class) for class in classes]

    #=check wonderful orthogonality theorem for character table=#
    for i in 1:irrep_num #eachindex(character_table[1])
        for j in 1:irrep_num #eachindex(character_table[1])
            result = sum([character_table[i, idx]*conj(character_table[j, idx])*multi[idx] for idx in eachindex(classes)])
            if i == j
                if !all(norm.(round.(result/gorder-complex(1.0,0.0), digits=5)) .< tor) return false end
            else
                if !all(norm.(round.(result/gorder, digits=5)) .< tor) return false end
            end
        end
    end
    #=check second orthogonality relation for character table=#
    for i in eachindex(classes)
        for j in eachindex(classes)
            result = sum([character_table[idx, i]*conj(character_table[idx, j]) * multi[i] for idx in 1:irrep_num])
            if i == j
                if !all(norm.(round.(result/gorder-complex(1.0,0.0), digits=5)) .< tor) return false end
            else
                if !all(norm.(round.(result/gorder, digits=5)) .< tor) return false end
            end
        end
    end
    return true
end

"""
get_character_table_of_irreps(irreps)
"""
function get_character_table_of_irreps(irreps; class::Bool, tor=1e-3)
    character_table = []
    gorder = size(irreps, 2)
    if class
        classes = []
        for i in 1:gorder
            table_coli = [tr(irrep) for irrep in irreps[:, i]]
            for j in axes(character_table, 2)
                if character_table == []
                    character_table = [table_coli;;]
                    push!(classes, [i])
                elseif all(norm.(round.(table_coli-character_table[:, j], digits=5)) .< tor)
                    push!(classes[j], i)
                else
                    continue
                end
                @goto next_i
            end
            character_table = [character_table;; table_coli]
            push!(classes, [i])
            @label next_i
        end
        return character_table, classes
    else
        character_table = tr.(irreps)
        return character_table
    end
end

end #module GroupRepresentation
