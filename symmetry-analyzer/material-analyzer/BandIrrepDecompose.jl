function read_structure(material_dir::String)
    INFO = string(material_dir, "/INFO_JULIA")
    POSCAR = string(material_dir, "/POSCAR_JULIA")
    return sys.get_structure(INFO, POSCAR)
end

function get_symmetries(cell, material_dir::String)
    println("Getting symmetries...")
    flush(stdout)
    syms = sys.get_symmetries(cell)
    mats = syms["operations"]
    
    open("$material_dir/symmetries.txt", "w") do file
        for mat in mats
            println(file, mat)
        end
    end
    println("Done. Symmetries saved.")
    return mats
end 

function get_irrep_decomposition_at_k(cell, material_dir::String; k)
    mats = []
    open("$material_dir/symmetries.txt") do file
        for line in readlines(file)
            push!(mats, eval(Meta.parse(line)))
        end
    end
    sys.get_irrep_decomposition_at_k(cell, mats; k)
end 