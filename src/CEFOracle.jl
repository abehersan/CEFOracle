module CEFOracle


using DataFrames
using LinearAlgebra
using StaticArrays
using Statistics
using OffsetArrays


include("./single_ion.jl")
export single_ion, mag_ion


include("./units.jl")
export meV_per_K, mu0, muB, kB, NA, Rg


include("./utils.jl")
export effective_moment
export is_normalized, is_hermitian, is_unitary


include("./blm_utils.jl")
export blm_dframe, alm_dframe
export get_blm!, get_alm!


include("./powder_grid.jl")
export cart_coords, columns
export SOPHE_xyzw
export SOPHE_grid


include("./cef_matrix.jl")
export cef_hamiltonian
export cef_eigensystem
export spin_operators
export stevens_EO


include("./thermodynamical_quantities.jl")
export population_factor
export partition_function


include("./cef_magnetization.jl")
export cef_magneticmoment_crystal!, cef_magneticmoment_powder!


include("./cef_susceptibility.jl")
export cef_susceptibility_crystal!, cef_susceptibility_powder!


include("./cef_entropy.jl")
export cef_entropy!, cef_entropy_speclevels!


include("./cef_neutron_xsection.jl")
export cef_neutronxsection_crystal!, cef_neutronxsection_powder!
export TAS_resfunc, gaussian, lorentz


include("./cef_ops.jl")
export stevens_O


include("./cef_rot.jl")
export rotate_blm
export get_euler_angles
export ZYZ_rotmatrix


end
