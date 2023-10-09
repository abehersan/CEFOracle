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

    """ INS X-section """
    temps = [10.0, 50.0, 200.0]
    ins_plot = plot(xlabel="Energy [meV]", ylabel="I(Q, E) [arb. units]")
    for t in temps
        calc_grid = DataFrame(T=t, Q=2.55, EN=range(-12,12,700))
        cef_neutronxsection_powder!(ion, bfactors, calc_grid)
        @df calc_grid plot!(ins_plot, :EN, :I_CALC, label="T=$t K")
    end

    """ Magnetic moment """
    mag_plot = plot(xlabel="Magnetic Field [Tesla]", ylabel="Magnetic Moment [muB]")
    calc_grid = DataFrame(T=1.5, Bx=0.0, By=0.0, Bz=0.0:0.5:12)
    cef_magneticmoment_crystal!(ion, bfactors, calc_grid, units="ATOMIC")
    @df calc_grid plot!(mag_plot, :Bz, :M_CALC, label="B parallel z")
    calc_grid = DataFrame(T=1.5, Bx=0:0.5:12, By=0.0, Bz=0.0)
    cef_magneticmoment_crystal!(ion, bfactors, calc_grid, units="ATOMIC")
    @df calc_grid plot!(mag_plot, :Bx, :M_CALC, label="B parallel x")
    ylims!(mag_plot, 0, 3.5)


    """ Static susceptibility """
    chi_plot = plot(xlabel="Temperature [K]", ylabel="1/chi [mol/cm^3]")
    calc_grid_z = DataFrame(T=1.5:0.5:300, Bx=0.0, By=0.0, Bz=0.01)
    cef_susceptibility_crystal!(ion, bfactors, calc_grid_z, units="CGS")
    @df calc_grid_z plot!(chi_plot, :T, 1 ./ :CHI_CALC, label="B parallel z")
    calc_grid_x = DataFrame(T=1.5:0.5:300, Bx=0.01, By=0.0, Bz=0.0)
    cef_susceptibility_crystal!(ion, bfactors, calc_grid_x, units="CGS")
    @df calc_grid_x plot!(chi_plot, :T, 1 ./ :CHI_CALC, label="B parallel x")
    chi_powd = 1/3*calc_grid_z.CHI_CALC + 2/3*calc_grid_x.CHI_CALC
    plot!(chi_plot, calc_grid_z.T, 1 ./ chi_powd, label="Powder")

    """ Specific heat capacity and magnetic entropy """
    s_plot = plot(xlabel="Temperature [K]", ylabel="S, Cp [J/mol/K]")
    calc_grid = DataFrame(T=0.5:0.5:300)
    cef_entropy!(ion, bfactors, calc_grid)
    @df calc_grid plot!(s_plot, :T, :CP_CALC, label="Cp")
    @df calc_grid plot!(s_plot, :T, :S_CALC, label="S")

    plot(ins_plot, mag_plot, chi_plot, s_plot, layout=(2, 2))
end


yb_hexagonalcef()
