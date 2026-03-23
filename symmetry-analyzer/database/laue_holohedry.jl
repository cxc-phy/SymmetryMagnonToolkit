

module LaueHolohedry

export LAUE1
export LAUE2M
export LAUEMMM
export LAUE4M
export LAUE4MMM
export LAUE3
export LAUE3M
export LAUE6M
export LAUE6MMM
export LAUEM3
export LAUEM3M
export TRICLI
export MONOCLI
export ORTHO
export TETRA
export TRIGO
export HEXA
export CUBIC

abstract type LAUE end
abstract type LAUE_NONE <: LAUE end
abstract type LAUE1 <: LAUE end
abstract type LAUE2M <: LAUE end
abstract type LAUEMMM <: LAUE end
abstract type LAUE4M <: LAUE end
abstract type LAUE4MMM <: LAUE end
abstract type LAUE3 <: LAUE end
abstract type LAUE3M <: LAUE end
abstract type LAUE6M <: LAUE end
abstract type LAUE6MMM <: LAUE end
abstract type LAUEM3 <: LAUE end
abstract type LAUEM3M <: LAUE end

abstract type HOLOHEDRY end
abstract type HOLOHEDRY_NONE <: HOLOHEDRY end
abstract type TRICLI <: HOLOHEDRY end
abstract type MONOCLI <: HOLOHEDRY end
abstract type ORTHO <: HOLOHEDRY end
abstract type TETRA <: HOLOHEDRY end
abstract type TRIGO <: HOLOHEDRY end
abstract type HEXA <: HOLOHEDRY end
abstract type CUBIC <: HOLOHEDRY end
end
