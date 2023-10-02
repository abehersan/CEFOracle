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
export blm_dframe, alm_dframe
export get_blm!, get_alm!


include("./cef_matrix.jl")
export cef_hamiltonian, cef_wavefunctions, cef_energies
export cef_eigensystem, spin_operators
# export cef_eigensystem_multisite, cef_site


include("./mag_properties.jl")
export cef_magnetization_crystal, cef_magnetization_powder #, cef_magnetization_multisite
export cef_susceptibility_crystal, cef_susceptibility_powder #, cef_susceptibility_multisite
export cef_heatcapacity, cef_heatcapacity_speclevels
export cef_entropy, cef_entropy_speclevels


include("./neutron_xsection.jl")
export cef_neutronxsection_crystal, cef_neutronxsection_powder #, cef_neutronxsection_multisite
export TAS_resfunc, gaussian, lorentz


include("./cef_ops.jl")
export stevens_O


include("./cef_rot.jl")
export rotate_blm, get_euler_angles, ZYZ_rotmatrix

# include("./cef_fit.jl")
# export cef_datasets, chi2_cef

end
