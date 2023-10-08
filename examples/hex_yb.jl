"""
Reproduction of results from the paper by Rotter and Bauer
"""


using CEFOracle
using DataFrames
using StatsPlots


ion = single_ion("Yb3")
bfactors = blm_dframe(Dict("B20"=>0.5622, "B40"=>1.6087e-5,
                           "B60"=>6.412e-7, "B66"=>-8.324e-6))


function yb_hexagonalcef()
    """
    calculate the inelastic neutron scattering x-section at different temperatures
    with model A, reproduces figure 7 of the paper
    """
    temps = [10.0, 50.0, 200.0]
    ins_plot = plot(xlabel="Energy [meV]", ylabel="I(Q, E) [arb. units]")
    for t in temps
        calc_grid = DataFrame(T=t, Q=2.55, EN=range(-12,12,700))
        cef_neutronxsection_powder!(ion, bfactors, calc_grid)
        @df calc_grid plot!(ins_plot, :EN, :CALC, label="T=$t K")
    end

    calc_grid = DataFrame(T=1.0:0.5:300)
    cef_entropy!()

    ins_plot
end


yb_hexagonalcef()