"""
units.jl

Define the `atomic`, `cgs` and `SI` system of units
"""


const meV_per_K = 0.086173332621451774  # divide En/meV_per_K to get En in K
# const mu0 = 201.33545383470705041      # T^2 Å^3 / meV
const mu0 = 1.25663706212e-6            # N/A
const muB = 0.057883818060738013331     # meV / T
const kB = 0.08617333262                # meV / K
const NA = 6.02214076e23                # Avogadro number mol^-1
const R = 8.314462618                   # Ideal gas constant J/mol/K


"""
TODO:
FINISH THIS SHITTTTTTTTTT
AT THE MOMENT 17.03, ALL CALCULATIONS OF PHYSICAL QUANTITIES ARE DONE
IN "SI UNITS" -> EVERYTHING HAS UNITS OF meV
"""
# Base.@kwdef struct physical_constants
#     mu0::Float64
#     muB::Float64
#     kB::Float64
#     NA::Float64
# end


# const units = (;
#     SI = physical_constants(;
#         mu0 = 201.33545383470705041,     # T^2 Å^3 / meV
#         muB = 0.057883818060738013331,   # meV / T
#         kB = 0.08617333262,             # meV / K
#         NA = 6.02214076e23              # Avogadro number mol^-1
#     ),
#     atomic = physical_constants(;
#         mu0 = 1.0,
#         muB = 1.0,
#         kB = 0.08617333262,
#         NA = 6.02214076e23
#     ),
#     cgs = physical_constants(;
#         mu0 = 1.0,
#         muB = 9.2740100783e-21,          # erg/G
#         kB = 0.08617333262,
#         NA = 6.02214076e23
#     ),
# )