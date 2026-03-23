function init()
    println("Initializing symmetry-analyzer and magnon-solver (~10 seconds)...")
    flush(stdout)
    cwd = pwd() #pathway for notebook.ipynb
    Base.include(Main, "$cwd/../symmetry-analyzer/material-analyzer/system.jl")
    Base.include(Main, "$cwd/../symmetry-analyzer/material-analyzer/SymmetryNeighborInfoPrep.jl")
    Base.include(Main, "$cwd/../symmetry-analyzer/material-analyzer/BandIrrepDecompose.jl")
    Base.include(Main, "$cwd/../magnon-solver/julia_fortran_interface.jl")
    println("Done.")
end