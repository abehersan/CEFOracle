"""
cef_multisite.jl

Define a multisite CEF matrix, diagonalize it and use the resulting
eigenvalues and eigenvectors to calculate physical observables
"""


using CEFOracle, Plots


ion = single_ion("Ho3")

bfactors1 = blm_dframe(Dict("B20"=>-0.089e-1))
bfactors2 = blm_dframe(Dict("B20"=>2.089e-3))

cef_site1 = cef_site(ion, bfactors1, 20.0, [1, 0, 0])
cef_site2 = cef_site(ion, bfactors2, 20.0, [0, 0, 1])

cef_sites = [cef_site1, cef_site2]
display(cef_eigensystem_multisite(cef_sites, verbose=true))

# temps = 1.9:0.01:100
# mag_powder = [cef_magnetization_multisite(cef_sites, T=t, units="ATOMIC") for t in temps]
# mag_p = plot(xlabel="Temperature [K]", ylabel="Calc. Magnetic Moment [B.M.]")
# plot!(mag_p, temps, mag_powder, label="Magnetization Ho3+ Multisite")

energy_transfer = -2:0.01:2
insx_powder = [cef_neutronxsection_multisite(cef_sites, en, 2.55, 2.0) for en in energy_transfer]
ins_p = plot(xlabel="Energy [meV]", ylabel="Calc. Intensity [arb. units]")
plot!(ins_p, energy_transfer, insx_powder, label="INS-XSECTION Ho3+ Multisite")

# display(mag_p)
display(ins_p)
