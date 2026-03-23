

module GroupIrrepresentationConstruct

include("kspace.jl")
include("grouprepresentation.jl")
include("../math-tool/mathtool.jl")
include("groupwithmultiplicationtable.jl")
include("../symmetry-group-tool/abstractgroup.jl")
include("../symmetry-group-tool/pointgroup.jl")
include("../symmetry-group-tool/spinpointgroup.jl")
include("../symmetry-operation-tool/symmetry.jl")

import .PointGroup as pg
import .SpinPointGroup as spg
import .SpaceSymmetry as spa
import .Symmetry
import .kSpace as kspa
import .GroupRepresentation as grep
import .GroupWithMultiplicationTable as gwmt
import .AbstractGroup as abg
using LinearAlgebra
using Combinatorics

function get_coirreps_of_little_group_of_k(mats; k)

    #=
    #= finite spin group=#
    Umats = [mat for mat in mats if det(mat.s) > 0]
    aUmats = [mat for mat in mats if det(mat.s) < 0]
    Umats, irreps = get_irreps_of_little_group_of_k(Umats, k=k);
    if length(aUmats) != 0
        totmats = [Umats; aUmats]
        coirrep = get_coirreps(totmats, Umats, aUmats, irreps, k=k, coirrep=true);
        uidxs = [i for i in eachindex(totmats) if Symmetry.is_symmetry_unitary(totmats[i])]
        irreps = [coirrep.irreps[i, j] for i in axes(coirrep.irreps, 1), j in uidxs];
    end
    =#

    #=
    #=enhanced spin group can find isomorphic abstract group=#
    mats = magspasym[change_symmetry_from_spin_basis_to_magnon_basis(mat, SO2=[0,0,1]) for mat in mats]
    mats = get_enhanced_group(mats)
    Umats = [mat for mat in mats if mat.T == 1]
    aUmats = [mat for mat in mats if mat.T == -1]
    Umats, irreps = get_irreps_of_little_group_of_k(Umats, k=k);
    if length(aUmats) != 0
        totmats = [Umats; aUmats]
        coirrep = get_coirreps(totmats, Umats, aUmats, irreps, k=k, coirrep=true);
        uidxs = [i for i in eachindex(totmats) if is_symmetry_unitary(totmats[i])]
        irreps = [coirrep.irreps[i, j] for i in axes(coirrep.irreps, 1), j in uidxs];
    end
    =#

    #=enhanced spin group can not find isomorphic abstract group=#
    mats = Symmetry.magspasym[Symmetry.change_symmetry_from_spin_basis_to_magnon_basis(mat, SO2=[0,0,1]) for mat in mats]
    Umats = [mat for mat in mats if mat.T == 1]
    aUmats = [mat for mat in mats if mat.T == -1]
    Pmats = [mat for mat in Umats if !all(abs.(diag(mat.s)) .< 1e-5)]#P'
    Qmats = [mat for mat in Umats if all(abs.(diag(mat.s)) .< 1e-5)]#Q'
    Pmats, irreps = get_irreps_of_little_group_of_k(Pmats, k=k)#irreps of P'
    Pmats, irreps = get_irreps_of_direct_product(Pmats, irreps)#irreps of Pgroup
    Umats = get_enhanced_group(Umats)#Ugroup
    aUmats = get_enhanced_group(aUmats)#aUgroup
    irreps = get_induced_irreps(Umats, Pmats, irreps, k=k)#induced irreps
    totmats = [Umats; aUmats]
    coirrep = get_coirreps(totmats, Umats, aUmats, irreps, k=k, coirrep=true);
    uidxs = [i for i in eachindex(totmats) if Symmetry.is_symmetry_unitary(totmats[i])]
    irreps = [coirrep.irreps[i, j] for i in axes(coirrep.irreps, 1), j in uidxs];
    return Umats, irreps
end

"""
REMARK!
    1.syms should be unitary nontrivial spin space group
    2.instead of true irreps (Γirreps), projective irreps (Dirreps) will be returned.
formula
    little group of k
        Hk = H1*T ⊕ H2*T ⊕ ... 
        T is the maximal translation group of Hk
        H = [H1, H2, ...] is the coset_representatives
    little co-group of k
        R = [(H1.r, H1.s), (H2.r, H2.s), ...] in spin basis
        or
        R = [(H1.r, H1.s, H1.T), (H2.r, H2.s, H2.T), ...] in collinear AFM magnon basis
    true irreps and projective irreps
        Γ(Hi) = D(Ri)*exp(-2π*im*k*Hi.t)
procedure
    1.get coset representatives: H
    2.get little co-group: R
    3.(1)if Gk is symmorphic or k inside BZ, get irrep D of point group R
    3.(2)otherwise, get projective irrep D of point group R, factor system decided by H
"""
function get_irreps_of_little_group_of_k(symmetries::Union{Vector{Symmetry.spinspasym}, Vector{Symmetry.magspasym}}; k, tor=1e-5)
    #=check if k is invariant under group=#
    syms = deepcopy(symmetries)
    cosetrep = Symmetry.get_coset_representatives_wrt_translation_group(syms)
    Dirreps = get_project_irrep_by_central_extension_method(cosetrep, k=k)
    return cosetrep, Dirreps
end

"""
module space group
get_project_irrep_by_central_extension_method(Hsyms::Vector{sg.sg_symtype}; k, tor=1e-5) -> Matrix{Matrix{ComplexF64}}
reference
    The Mathematical Theory of Symmetry in Solids: Representation Theory for Point Groups and Space Groups
    (C. J. Bradley, A. P. Cracknell)
    p155-p161
formula
    group product of project irrepresentation D(Ri) of R, decided by Hsyms
        D(Ri)*D(Rj) = μ(Ri, Rj)*D(RiRj)
        μ(Ri, Rj) = exp(-2π*im*Glist[i]*Hj.t)
        Glist[i] = Hi^{-1}*k - k
        Hi.r = Ri

    multiplier of H: Zh
        Zh = [0, 1, ..., h-1] where h is the minimal integer s.t. μ(Hi, Hj)^h = I
    group product
        i*j = (i+j) mod h

    central_extension of H with kernel Zh: H*
        H* = [(Hi, j)...] where Hi ∈ H and j ∈ Zh
    group product
        (Hi, j)*(Hk, l) = (HiHk, j+l+a(HiHk))
        a(Hi, Hj) = -Glist[i]*Hj.t*g -> μ(Hi, Hj) = exp(2π*im*a(Hi,Hj))
procedure
    1.get multiplier of H: Zh
    2.get central_extension of H with kernel Zh: G*
    3.get abstract group isomorphic to G*
    4.get irreps
"""
function get_project_irrep_by_central_extension_method(Hsyms::Union{Vector{Symmetry.spinspasym}, Vector{Symmetry.magspasym}}; k, tor=1e-5)#::Matrix{Matrix{ComplexF64}}
    #=get G = Hsym^-1*k -k for each Hsym =#
    Glist = []
    Ididx = 0
    for i in eachindex(Hsyms)
        sym = Hsyms[i]
        if Ididx == 0 && Symmetry.is_Identity(sym) Ididx = i end
        ksym = Symmetry.get_symmetry_in_kspace(sym)
        G = round.(inv(ksym)*k-k,digits=5)
        push!(Glist, G[1:3, 4])
    end
    #=get factor μsystem=#
    μsys = zeros(ComplexF64, length(Hsyms), length(Hsyms))
    for i in eachindex(Hsyms), j in eachindex(Hsyms)
        μsys[i, j] = exp(-2.0*π*im*sum(Glist[i].*Hsyms[j].t))
    end
    #=get number h=#
    h = 0 
    for i in 1:10
        if !all(norm.(round.(μsys .^i.-1, digits=2)) .< tor) continue end
        h = i
        break
    end
    if h == 0 error("not find number h") end
    #=get asystem=#
    asys = zeros(Int64, length(Hsyms), length(Hsyms))
    for i in eachindex(Hsyms), j in eachindex(Hsyms)
        asys[i, j] = Int64.(round.(log(μsys[i, j])*h/(2.0*π*im), digits=3))
    end
    #=get multiplication table of central extension=#
    if typeof(Hsyms) == Vector{Symmetry.spinspasym}
        syms = [(r=sym.r, s=sym.s) for sym in Hsyms]
    else
        syms = [(r=sym.r, s=sym.s, T=sym.T) for sym in Hsyms]
    end
    multable = Symmetry.get_multiplication_table_of_group(syms)
    ext_syms = [(i, j) for i in eachindex(syms) for j in 0:h-1]
    ext_multable = get_multiplication_table_of_central_extension(ext_syms, multable, h, asys)
    classes = gwmt.get_full_conjugacy_classes(multable=ext_multable)
    #=get irrep of central extension=#
    idx, perm = get_class_of_abstract_group_and_isomorphism(ext_multable, classes)
    #println(idx)
    #println(abg.abg_data[idx].notation)
    irrep = get_irreps_of_abstract_group(idx)
    ext_irrep = [ irrep[i, perm[j]] for i in axes(irrep, 1), j in axes(irrep, 2)]
    #=get allowed irreps of central extension=#
    allowed_irrep = []
    for i in axes(ext_irrep, 1)
        for j in 0:h-1
            if !all(norm.(round.(ext_irrep[i, h*(Ididx-1)+j+1] - exp(2.0*π*im*j/h)*I, digits=5)) .< tor) @goto next_i end
        end
        if allowed_irrep == []
            allowed_irrep = reshape(ext_irrep[i,:], 1, size(ext_irrep, 2))
        else
            allowed_irrep = [allowed_irrep;reshape(ext_irrep[i,:], 1, size(ext_irrep, 2))]
        end
        @label next_i     
    end
    #=get projective irreps of R=#
    proirreps = []
    for i in axes(allowed_irrep, 2)
        if ext_syms[i][2] != 0 continue end
        if proirreps == []
            proirreps = [allowed_irrep[:, i];;]
        else
            proirreps = [proirreps ;; allowed_irrep[:, i]]
        end
    end
    return proirreps
end

""" 
module space group
get_multiplication_table_of_central_extension(ext_syms, multable, h, asys) -> Matrix{Int64}
"""
function get_multiplication_table_of_central_extension(ext_syms::Vector{Tuple{Int64, Int64}}, multable::Matrix{Int64}, h::Int64, asys::Matrix{Int64})::Matrix{Int64}
    order = length(ext_syms)
    ext_multable = zeros(Int64, order, order)
    for i = 1:order
        for j = 1:order
            sym1 = ext_syms[i]
            sym2 = ext_syms[j]
            nsym = (multable[sym1[1],sym2[1]], mod(sym1[2]+sym2[2]+asys[sym1[1],sym2[1]], h))
            idx = findall(sym->all(sym.==nsym), ext_syms)[1]
            ext_multable[i, j] = idx
        end
    end
    return ext_multable
end

"""
module abstract group
get_class_of_abstract_group_and_isomorphism(multable::Matrix{Int64}, classes::Vector{Vector{Int64}})
"""
function get_class_of_abstract_group_and_isomorphism(multable::Matrix{Int64}, classes::Vector{Vector{Int64}})
    genset = gwmt.get_a_minimal_generating_set_of_group(multable=multable)
    orderlist = [gwmt.get_order_of_symmetry(gen, multable=multable) for gen in genset]
    syms = Vector{Int64}[]
    for idx in axes(multable, 1)
        sym = gwmt.sort_symmetry_in_order_of_generators(idx, gens=genset, orderlist=orderlist, multable=multable)
        if sym == [] error(0) end
        push!(syms, sym)
    end
    inf = sort([length(class) for class in classes])
    for i in eachindex(abg.abg_multiplication_table)
        itable = abg.abg_multiplication_table[i]
        if itable == [] continue end
        if length(syms) != size(itable, 1) continue end
        iclasses = gwmt.get_full_conjugacy_classes(multable=itable)
        iinf = sort([length(class) for class in iclasses])
        if inf != iinf continue end
        permutation = get_isomorphism(multable, genset, orderlist, syms, itable)
        if permutation != [] return i, permutation end
    end
    error("can not find class of abstract group")
end

"""
get_isomorphism(multable1::Matrix{Int64}, gens1, genorderlist1, syms1, multable2::Matrix{Int64}) -> Vector{Int64}
    permutation:P[i]
    syms1 = [1,2,3,...,i]
    syms2 = [,....,P[i],]
REMARK!
    
"""
function get_isomorphism(multable1::Matrix{Int64}, gens1::Vector{Int64}, genorderlist1::Vector{Int64}, syms1::Vector{Vector{Int64}}, multable2::Matrix{Int64})::Vector{Int64}
    #=
    if size(multable1, 1) != size(multable2, 1) return [] end
    classes1 = gwmt.get_full_conjugacy_classes(multable=multable1)
    inf1 = sort([length(class) for class in classes1])
    classes2 = gwmt.get_full_conjugacy_classes(multable=multable2)
    inf2 = sort([length(class) for class in classes2])
    if inf1!=inf2 return [] end
    =#
    
    #=get minimal generating sets of multable2 and compare=#
    len = length(gens1)
    orderlist1 = deepcopy(genorderlist1)
    Eidx = gwmt.get_index_of_identity_element(multable=multable2)
    idxs = [i for i in axes(multable2, 1) if i != Eidx]
    if length(idxs) == 0 return [1] end
    for perm in permutations(idxs, len)
        if length(Set(gwmt.generate_full_elements_in_order_of_generators(perm, multable=multable2))) != size(multable2, 1)
            continue
        end
        orderlist2 = [gwmt.get_order_of_symmetry(gen, multable=multable2) for gen in perm]
        if orderlist2 != orderlist1 continue end
        syms2 = Vector{Int64}[]
        for idx in axes(multable2, 1)
            sym2 = gwmt.sort_symmetry_in_order_of_generators(idx, gens=perm, orderlist=orderlist2, multable=multable2)
            if sym2 == [] error(2) end
            push!(syms2, sym2)
        end
        if !gwmt.is_group_same(syms1, syms2) continue end
        permutation = Vector{Int64}([findall(x->x==sym, syms2)[1] for sym in syms1])
        if all(multable2-gwmt.permute_multiplication_table(permutation, multable1).==0)
            return permutation
        end
    end
    return Int64[]
end

"""
module abstract group
get irreps of abstract group
"""
function get_irreps_of_abstract_group(idx::Int64)::Matrix{Matrix{ComplexF64}}
    irreps = abg.abg_irreps[idx]
    row = length(irreps)
    col = length(irreps[1])
    irreps = [Matrix{ComplexF64}([irreps[i][j];;]) for i in 1:row, j in 1:col]
    if !grep.check_irreps(irreps) error("wrong irreps of abstract group") end
    return irreps
end
function get_enhanced_group(syms::Vector{Symmetry.spinspasym};SO2)::Vector{Symmetry.magspasym}
    G = Symmetry.magspasym[
        (r=[1 0 0; 0 1 0; 0 0 1], t=[0.0, 0.0, 0.0], s=[ 1.0 0.0; 0.0  1.0], T=1),
        (r=[1 0 0; 0 1 0; 0 0 1], t=[0.0, 0.0, 0.0], s=[ 1.0 0.0; 0.0 -1.0], T=1),
        (r=[1 0 0; 0 1 0; 0 0 1], t=[0.0, 0.0, 0.0], s=[-1.0 0.0; 0.0 -1.0], T=1),
        (r=[1 0 0; 0 1 0; 0 0 1], t=[0.0, 0.0, 0.0], s=[-1.0 0.0; 0.0  1.0], T=1)
    ]
    nsyms = Symmetry.magspasym[]
    for sym in syms
        #println(1,sym)
        nsym = Symmetry.change_symmetry_from_spin_basis_to_magnon_basis(sym, SO2=SO2)
        #println(2,nsym)
        for g in G
            push!(nsyms, Symmetry.composite_symmetries(g, nsym, modtrans=true))
        end
    end
    if !Symmetry.isa_group(nsyms) error("not a group") end
    return nsyms
end
function get_enhanced_group(syms::Vector{Symmetry.magspasym})::Vector{Symmetry.magspasym}
    G = Symmetry.magspasym[
        (r=[1 0 0; 0 1 0; 0 0 1], t=[0.0, 0.0, 0.0], s=[ 1.0 0.0; 0.0  1.0], T=1),
        (r=[1 0 0; 0 1 0; 0 0 1], t=[0.0, 0.0, 0.0], s=[ 1.0 0.0; 0.0 -1.0], T=1),
        (r=[1 0 0; 0 1 0; 0 0 1], t=[0.0, 0.0, 0.0], s=[-1.0 0.0; 0.0 -1.0], T=1),
        (r=[1 0 0; 0 1 0; 0 0 1], t=[0.0, 0.0, 0.0], s=[-1.0 0.0; 0.0  1.0], T=1)
    ]
    nsyms = Symmetry.magspasym[]
    for sym in syms
        for g in G
            push!(nsyms, Symmetry.composite_symmetries(g, sym, modtrans=false))
        end
    end
    #if !Symmetry.isa_group(nsyms) error("not a group") end
    return nsyms
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
function get_irreps_of_direct_product(syms::Vector{Symmetry.magspasym}, irreps::Matrix{Matrix{Complex{Float64}}})#::Matrix{Matrix{Complex{Float64}}}
    G = Symmetry.magspasym[
        (r=[1 0 0; 0 1 0; 0 0 1], t=[0.0, 0.0, 0.0], s=[ 1.0 0.0; 0.0  1.0], T=1),
        (r=[1 0 0; 0 1 0; 0 0 1], t=[0.0, 0.0, 0.0], s=[ 1.0 0.0; 0.0 -1.0], T=1),
        (r=[1 0 0; 0 1 0; 0 0 1], t=[0.0, 0.0, 0.0], s=[-1.0 0.0; 0.0 -1.0], T=1),
        (r=[1 0 0; 0 1 0; 0 0 1], t=[0.0, 0.0, 0.0], s=[-1.0 0.0; 0.0  1.0], T=1)
    ]
    Girreps = [
        #1  1  1  1
        #1 -1  1 -1
        1  1 -1 -1
        1 -1 -1  1
    ]
    
    #=get irreps of direct product of nspg and {[E||E], [Σz||E], [-E||E], [-Σz||E]}=#
    nsyms = Symmetry.magspasym[]
    for g in G, sym in syms
        push!(nsyms, Symmetry.composite_symmetries(g, sym, modtrans=false))
    end
    directproduct_irreps = Matrix{Any}([nothing for i=1:2*size(irreps, 1), j=1:4*size(irreps, 2)])
    for i in axes(Girreps, 1), j in axes(Girreps, 2)
        directproduct_irreps[(i-1)*size(irreps, 1)+1:i*size(irreps, 1), (j-1)*size(irreps, 2)+1:j*size(irreps, 2)] = Girreps[i,j]*irreps
    end
    
    
    for i in axes(directproduct_irreps, 1), j in axes(directproduct_irreps, 2)
        if size(directproduct_irreps[i, j]) == () directproduct_irreps[i, j] = [directproduct_irreps[i, j];;] end
    end
    directproduct_irreps = Matrix{Matrix{Complex{Float64}}}(directproduct_irreps)
    return nsyms, directproduct_irreps
end

"""
get induced irrep of super from sub(only solve invariant subgroup of index 2)
input
    D(Ri)
output
    D(Ri)
true irreps and projective irreps
    Γ(Hi) = D(Ri)*exp(-2π*im*k*Hi.t)
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

function get_induced_irreps(super::Vector{Symmetry.magspasym}, sub::Vector{Symmetry.magspasym}, irreps_list::Matrix{Matrix{Complex{Float64}}}; k, tor=1e-3)#::Matrix{Matrix{Complex{Float64}}}
    if length(super) == length(sub) 
        idxs = [findall(x->Symmetry.is_symmetry_same(x, mat, modtrans=false), sub)[1] for mat in super]
        irreps_list = [irreps_list[i, j] for i in axes(irreps_list, 1), j in idxs]
        return irreps_list 
    end
    if length(super) != 2* length(sub) error("only solve invariant subgroup of index 2") end

    q = [mat for mat in super if !Symmetry.is_symmetry_contained(mat, sub, modtrans=false)][1]#coset representative
    conj_sub = [Symmetry.conjugate_symmetry(mat, q, modtrans=false) for mat in sub]
    conj_idx = [findall(x->Symmetry.is_symmetry_same(mat, x, modtrans=true), sub)[1] for mat in conj_sub]####mention modtrans=true
    gorder = size(irreps_list, 2)
    used_irrep_idx = []
    ind_irreps_list = []
    irreps_list = [irreps_list[i, idx]*exp(-2.0*π*im*sum(k[:,4].*sub[idx].t)) for i in axes(irreps_list, 1), idx in eachindex(sub)]
    for i in axes(irreps_list, 1)
        if i in used_irrep_idx continue end
        #=get conjugate irrepresentation=#
        conj_irreps = [irreps_list[i, idx] for idx in conj_idx]
        #=induce irrep depending on if irrep and conj_irrep are similar=#
        j = [m for m in axes(irreps_list, 1) if all(abs.(tr.(conj_irreps) - tr.(irreps_list[m, :])) .< tor)][1]
        dim = size(irreps_list[i, 1], 1)
        #println(i,j)
        if i == j
            Tmatrix = grep.get_similar_transformation_matrix_between_two_reps(irreps_list[i, :], conj_irreps)
            q_irreps = [Tmatrix, -Tmatrix]
            ind_irreps = Matrix{Matrix{Complex{Float64}}}([zeros(dim, dim) for i in 1:2, j in 1:length(super)])
            for pmidx = 1:2
                q_irrep = q_irreps[pmidx]
                for m in eachindex(sub)
                    u = sub[m]
                    u_irrep = Matrix{Complex{Float64}}(irreps_list[i, m])
                    Idx = findall(x->Symmetry.is_symmetry_same(x, u, modtrans=false), super)[1]
                    ind_irreps[pmidx, Idx] = u_irrep
                    qu = Symmetry.composite_symmetries(q, u, modtrans=false)
                    qu_irrep = MathTool.composite_matrice(q_irrep, u_irrep)
                    Idx = findall(x->Symmetry.is_symmetry_same(x, qu, modtrans=true), super)[1]
                    ind_irreps[pmidx, Idx] = qu_irrep*exp(-2.0*π*im*sum(k[:,4].*(qu.t-super[Idx].t)))
                end
            end
        else
            # not modified
            q2 = Symmetry.composite_symmetries(q, order=2, modtrans=false)
            q2_idx = findall(x->Symmetry.is_symmetry_same(q2, x, modtrans=true), sub)[1]
            q_irrep = Matrix{Complex{Float64}}([[zeros(dim, dim); Matrix(I, dim, dim)];; [irreps_list[i, q2_idx]*exp(-2.0*π*im*sum(k[:,4].*(q2.t-super[q2_idx].t))); zeros(dim, dim)]])
            ind_irreps = Matrix{Matrix{Complex{Float64}}}([zeros(2*dim, 2*dim) for i in 1:1, j in 1:length(super)])
            for m in eachindex(sub)
                u = sub[m]
                u_irrep = Matrix{Complex{Float64}}([[irreps_list[i, m]; zeros(dim, dim)];; [zeros(dim, dim); conj_irreps[m]]])
                Idx = findall(x->Symmetry.is_symmetry_same(x, u, modtrans=false), super)[1]
                ind_irreps[1, Idx] = u_irrep
                qu = Symmetry.composite_symmetries(q, u, modtrans=false)
                qu_irrep = MathTool.composite_matrice(q_irrep, u_irrep)
                Idx = findall(x->Symmetry.is_symmetry_same(x, qu, modtrans=true), super)[1]
                ind_irreps[1, Idx] = qu_irrep*exp(-2.0*π*im*sum(k[:,4].*(qu.t-super[Idx].t)))
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
    ind_irreps_list = [ind_irreps_list[i,j]*exp(2.0*π*im*sum(k[:,4].*super[j].t)) for i in axes(ind_irreps_list, 1), j in axes(ind_irreps_list, 2)]
    return ind_irreps_list
end

"""
#=non-unitary group, get coirrepresentation=#
"""
function get_coirreps(fullgroup::Union{Vector{Symmetry.magspasym}, Vector{Symmetry.spinspasym}}, Ugroup::Union{Vector{Symmetry.magspasym}, Vector{Symmetry.spinspasym}}, aUgroup::Union{Vector{Symmetry.magspasym}, Vector{Symmetry.spinspasym}}, irreps::Matrix{Matrix{ComplexF64}}; k, coirrep::Bool)#::NamedTuple{(:irreps, :case), Tuple{Matrix{Matrix{ComplexF64}}, Vector{Int64}}}#Tuple{Matrix{Matrix{ComplexF64}}, Vector{Int64}}
    #=get Wigner criterion=#
    character_sum = get_case_of_irreps(irreps, Ugroup, aUgroup, k)    
    if !coirrep
        error(111)
    else
        #=get coirreps=#
        case = Int64[]
        coirreps = Matrix{ComplexF64}[]
        a0 = aUgroup[1]
        a0_inv = Symmetry.get_inverse_of_symmetry(a0, modtrans=true)
        used_idx = []
        for i in axes(irreps, 1)
            if i in used_idx continue end
            conj_irrep = []
            for u in Ugroup
                u_a0_conj = Symmetry.composite_symmetries(a0_inv, u, a0, modtrans=false)
                u_a0_conj_idx = findall(x->Symmetry.is_symmetry_same(x, u_a0_conj, modtrans=true), Ugroup)[1]
                push!(conj_irrep, conj(irreps[i, u_a0_conj_idx]))
            end
            if size(irreps[i, 1]) == ()
                zeromat = ComplexF64(0)
            else
                zeromat = zeros(ComplexF64, size(irreps[i, 1]))
            end
            coirrepi = Matrix{ComplexF64}[]
            if character_sum[i] == 1
                a0mat = grep.get_similar_transformation_matrix_between_two_reps(irreps[i,:], conj_irrep)
                for sym in fullgroup
                    if Symmetry.is_symmetry_unitary(sym)
                        idx = findall(x->Symmetry.is_symmetry_same(x, sym, modtrans=false), Ugroup)[1]
                        push!(coirrepi, irreps[i, idx])
                    else
                        a = sym
                        aa0_inv = Symmetry.composite_symmetries(a, a0_inv, modtrans=false)
                        idx = findall(x->Symmetry.is_symmetry_same(x, aa0_inv, modtrans=true), Ugroup)[1]
                        amat = irreps[i, idx] * a0mat
                        push!(coirrepi, amat)
                    end
                end
            elseif character_sum[i] == 0
                push!(used_idx, grep.get_index_of_irrep(conj_irrep, irreps))
                for sym in fullgroup
                    if Symmetry.is_symmetry_unitary(sym)
                        u = sym
                        idx1 = findall(x->Symmetry.is_symmetry_same(x, u, modtrans=false), Ugroup)[1]
                        u_a0_conj = Symmetry.composite_symmetries(a0_inv, u, a0, modtrans=false)
                        idx2 = findall(x->Symmetry.is_symmetry_same(x, u_a0_conj, modtrans=true), Ugroup)[1]
                        umat = [[irreps[i, idx1] ;; zeromat] ; [zeromat ;; conj(irreps[i, idx2])]]
                        push!(coirrepi, umat)
                    else
                        a = sym
                        aa0 = Symmetry.composite_symmetries(a, a0, modtrans=false)
                        idx1 = findall(x->Symmetry.is_symmetry_same(x, aa0, modtrans=true), Ugroup)[1]
                        a0_inv_a = Symmetry.composite_symmetries(a0_inv, a, modtrans=false)
                        idx2 = findall(x->Symmetry.is_symmetry_same(x, a0_inv_a, modtrans=true), Ugroup)[1]
                        amat = [[zeromat ;; irreps[i, idx1]] ; [conj(irreps[i, idx2]) ;; zeromat]]
                        push!(coirrepi, amat)
                    end
                end
            elseif character_sum[i] == -1
                a0mat = grep.get_similar_transformation_matrix_between_two_reps(irreps[i,:], conj_irrep)
                for sym in fullgroup
                    if Symmetry.is_symmetry_unitary(sym)
                        idx = findall(x->Symmetry.is_symmetry_same(x, sym, modtrans=false), Ugroup)[1]
                        umat = [[irreps[i, idx] ;; zeromat] ; [zeromat ;; irreps[i, idx]]]
                        push!(coirrepi, umat)
                    else
                        a = sym
                        a_a0_inv = Symmetry.composite_symmetries(a, a0_inv, modtrans=false)
                        idx = findall(x->Symmetry.is_symmetry_same(x, a_a0_in, modtrans=true), Ugroup)[1]
                        mat = irreps[i, idx] * a0mat
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
        #coirreps = [coirreps[i, j]*exp(-2.0*π*im*sum(k[:,4].*fullgroup[j].t)) for i in axes(coirreps, 1), j in axes(coirreps, 2)]
        return (irreps=coirreps, case=case)
    end
end

function get_case_of_irreps(irreps, Ugroup, aUgroup, k)
    character_sum = zeros(ComplexF64, size(irreps, 1))
    for sym in aUgroup
        sym_squared = Symmetry.composite_symmetries(sym, order=2, modtrans=false)
        trans = sym_squared.t
        idx = findall(x->Symmetry.is_symmetry_same(x, sym_squared, modtrans=true), Ugroup)[1]
        for i = 1:size(irreps, 1)
            character_sum[i] += tr(irreps[i, idx])*exp(-2.0*π*im*sum(k[:,4].*trans))/length(Ugroup)
        end
    end
    character_sum = round.(character_sum, digits=3)
    if all(isinteger.(character_sum))
        character_sum = Int64.(character_sum)
    else
        println(character_sum)
        error("Error: character sum is not a real integer")
    end
    return character_sum
end

end #module GroupIrrepresentationConstruct