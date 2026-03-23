import .System as sys

function generate_neighobor_info(material_dir::String)
    cwd = pwd() #pathway for notebook.ipynb
    mat_path = joinpath(cwd, material_dir)
    cell = sys.get_structure("$mat_path/INFO_JULIA", "$mat_path/POSCAR_JULIA")
    infofort = sys.read_INFO_for_FORTRAN("$mat_path/INFO_FORT")
    nbr_order = infofort["nbr_order"]
    sym = infofort["sym"]
    if sym
        println("generate neighbor information to $(nbr_order)th order(symmetry: on)...")
    else
        println("generate neighbor information to $(nbr_order)th order(symmetry: off)...")
    end

    
    if cell["mag"]
        atom, idx_offset = sys.get_atom_and_index_offset_of_site_with_finite_magmom(cell)
        neighbors = sys.get_neighbors_between_atoms(cell, 1:nbr_order, atom, atom, sym=sym)
    else
        idx_offset = 0
        neighbors = sys.get_neighbors(cell, 1:nbr_order, sym=sym)
    end
    
    open("$mat_path/neighbor_info.txt", "w") do file
        #=first line of neighbor_info.txt=#
        if sym
            print(file, "sym: on, distances:")
        else
            print(file, "sym: off, distances:")
        end
        for i in sort(collect(keys(neighbors)))
            print(file, " ", neighbors[i]["dist"])
        end
        println(file)
        #=neighbor information=#
        for i in sort(collect(keys(neighbors)))
            ineq_neighbor = [ key for key in collect(keys(neighbors[i])) if typeof(key)==Int64]
            for j in ineq_neighbor#sort(collect(keys(neighbors[i])))
                for bond in neighbors[i][j]
                    print(file, i, " ", j, " ")
                    atom1 = bond[1]-idx_offset
                    atom2 = bond[2]-idx_offset
                    cell_dist = bond[3]
                    print(file, atom1, " ", atom2, " ", cell_dist[1], " ", cell_dist[2], " ", cell_dist[3], " ")
                    if cell["mag"]
                        println(file)
                    else
                        dir_cosine = MathTool.get_direction_cosine_of_vector(sys.get_cartesian_vector_between_sites(cell, atom1, atom2, cell_dist), x, y, z)
                        println(file, dir_cosine[1], " ", dir_cosine[2], " ", dir_cosine[3])
                    end
                end
            end
        end
    end
    println("done")
end


#=
cell = sys.get_structure("./full-INFO", "./full-POSCAR")
infofort = sys.read_INFO_for_FORTRAN("./INFO_FORT")
nbr_order = infofort["nbr_order"]
sym = infofort["sym"]

if cell["mag"]
    atom, idx_offset = sys.get_atom_and_index_offset_of_site_with_finite_magmom(cell)
    neighbors = sys.get_neighbors_between_atoms(cell, 1:nbr_order, atom, atom, sym=sym)
else
    idx_offset = 0
    neighbors = sys.get_neighbors(cell, 1:nbr_order, sym=sym)
end
#z=[0,0,1.0]
#x=cell["lattice"]*[1,1,0]
#y=cell["lattice"]*[-1,1,0]

open("neighbor_info.txt", "w") do file
    #=first line of neighbor_info.txt=#
    if sym
        print(file, "sym: on, distances:")
    else
        print(file, "sym: off, distances:")
    end
    for i in sort(collect(keys(neighbors)))
        print(file, " ", neighbors[i]["dist"])
    end
    println(file)
    #=neighbor information=#
    for i in sort(collect(keys(neighbors)))
        ineq_neighbor = [ key for key in collect(keys(neighbors[i])) if typeof(key)==Int64]
        for j in ineq_neighbor#sort(collect(keys(neighbors[i])))
            for bond in neighbors[i][j]
                print(file, i, " ", j, " ")
                atom1 = bond[1]-idx_offset
                atom2 = bond[2]-idx_offset
                cell_dist = bond[3]
                print(file, atom1, " ", atom2, " ", cell_dist[1], " ", cell_dist[2], " ", cell_dist[3], " ")
                if cell["mag"]
                    println(file)
                else
                    dir_cosine = MathTool.get_direction_cosine_of_vector(sys.get_cartesian_vector_between_sites(cell, atom1, atom2, cell_dist), x, y, z)
                    println(file, dir_cosine[1], " ", dir_cosine[2], " ", dir_cosine[3])
                end
            end
        end
    end
end
println("done")
=#