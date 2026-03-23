using DelimitedFiles
using Plots

function run_magnon_solver(material_dir::String)
    cwd = pwd() #pathway for notebook.ipynb
    exe_path = joinpath(cwd, "../magnon-solver/exe/magnon_solver.exe")
    mat_path = joinpath(cwd, material_dir)
    cd(mat_path)
    run(`$exe_path`)
    cd(cwd)
end

function plot_magnon_dispersion(material_dir::String; klabel=nothing)
    data = readdlm("$material_dir/dispersion.txt")
    num_bands = size(data, 2) - 4
    start_band = div(num_bands, 2) + 1
    k = Int.(data[:, 1])
    E = data[:, start_band+1:num_bands+1]

    tick_positions = data[1, 1]:5000:data[end, 1]+1
    if klabel == nothing
        tick_labels = fill("", length(tick_positions))
    else
        if length(tick_positions) != length(klabel)
            println("warning: number of labels doesn't match the number of high-symmetry points; should be $(length(tick_positions))." )
            tick_labels = fill("", length(tick_positions))
        else
            tick_labels = klabel
        end
    end
    plot(
        k, 
        E[:, 1],
        legend = false,
        xticks = (tick_positions, tick_labels),
        xrange = [data[1, 1],data[end, 1]],
        yrange = [
            maximum(
                [0, minimum(data[:,start_band+1])-0.1]
                ), 
            maximum(data[:,num_bands+1])+0.5],
        framestyle=:box
    )
    for i in 2:size(E, 2)
        plot!(k, E[:,i])
    end
    xlabel!("kpath")
    ylabel!("Energy")
end