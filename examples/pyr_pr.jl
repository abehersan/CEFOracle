"""
pyr_pr.jl

Calculate the inelastic neutron scattering cross-section with the parameters
given in Boothroyd's paper on Pr2Sn2O7
"""


using CEFOracle, LaTeXStrings, Plots
gr()


pr = single_ion("Pr3")
bs_dict = Dict("B20"=>57.9*1e-5,
               "B40"=>432.8*1e-5,
               "B43"=>161.0*1e-5,
               "B60"=>144.5*1e-5,
               "B63"=>-107.5*1e-5,
               "B66"=>192.3*1e-5)
bs_dframe = blm_dframe(bs_dict)

display(cef_eigensystem(pr, bs_dframe, verbose=true))

es = LinRange(-20, 100, 700)
temps = [5.0, 25.0, 100.0]
ins_plot = plot(xlabel="Energy [meV]", ylabel="I(Q, E) [arb. units]")
for t in temps
    plot!(es, [cef_neutronxsection(pr, bs_dframe, e, 2.55, t) for e in es],
          label="T=$t")
end
display(ins_plot)
# plot!(es, ins_xsection, label="T=5K")