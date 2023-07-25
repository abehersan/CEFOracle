"""
cef_multisite.jl

Define a multisite CEF matrix, diagonalize it and use the resulting
eigenvalues and eigenvectors to calculate physical observables
"""


using CEFOracle, Plots


ion = single_ion("Ho3")
bfactors1 = blm_dframe(Dict("B20"=>-0.089e-1))
bfactors2 = blm_dframe(Dict("B20"=>-2.089e-3))
sites_vecB = [cef_site(ion, bfactors1,0.5, [1,0,0]), cef_site(ion, bfactors2,0.5, [0,1,0])]
cef_eigensystem_multisite(sites_vecB, verbose=true);


energy_transfer = -2:0.01:2
insx_powder = [cef_neutronxsection_multisite(sites_vecB, en, 2.55, 2.0) for en in energy_transfer]


p = plot(xlabel="Energy [meV]", ylabel="Calc. Intensity [arb. units]")
plot(p, energy_transfer, insx_powder, label="INS-XSECTION Ho3+ Multisite")
