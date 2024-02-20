module CEFOracle


using DataFrames
using DataFramesMeta
using LinearAlgebra
using StaticArrays
using Statistics
using OffsetArrays


const PREC::Float64 = 1.0e-9    # for degeneracy calculations
const SDIG::Int64 = 9           # for numerical cutoffs


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
export thermal_average


include("./cef_magnetization.jl")
export cef_magneticmoment_crystal!
export cef_magneticmoment_powdergrid!


include("./cef_susceptibility.jl")
export cef_susceptibility_crystal!
export cef_susceptibility_powder!
export cef_susceptibility_powdergrid!


include("./cef_entropy.jl")
export cef_entropy!, cef_entropy_speclevels!


include("./cef_neutronxsection.jl")
export cef_neutronxsection_crystal!, cef_neutronxsection_powder!
export TAS_resfunc, gaussian, lorentz


include("./cef_ops.jl")
export stevens_O


include("./cef_rot.jl")
export rotate_blm
export get_euler_angles
export ZYZ_rotmatrix


end