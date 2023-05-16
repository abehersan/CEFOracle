module CEFOracle


# using external modules
using DataFrames
using JuMP
using LinearAlgebra
using Logging
using NLopt
using Plots
using SpecialFunctions
using StaticArrays
using Statistics
using OffsetArrays


# export functions and variables for REPL use
export single_ion
export meV_per_K, mu0, muB, kB, NA, R
export effective_moment, blm_dframe
export cef_eigensystem
export cef_magnetization, cef_susceptibility
export cef_heatcapacity, cef_heatcapacity_speclevels
export TAS_res, voigt, gauss, lorentz, cef_neutronxsection


# include source files in module scope
include("./single_ion.jl")
include("./units.jl")
include("./utils.jl")
include("./cef_matrix.jl")
include("./mag_properties.jl")
include("./neutron_xsection.jl")
include("./cef_fit.jl")
include("./cef_ops.jl")


end