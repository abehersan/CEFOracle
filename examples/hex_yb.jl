"""
Reproduction of results from Rotter, Bauer paper, CEF model A
"""


using CEFOracle
using GLMakie
using LaTeXStrings


yb = single_ion("Yb3")
bs = Dict("B20"=>0.5622, "B40"=>1.6087e-5, "B60"=>6.412e-7, "B66"=>-8.324e-6)
bdf = blm_dframe(bs)


"""
reproduces figure 11, page 220
"""
function plot_mag(ion, blm)
    T = 0.5
    fields = LinRange(0.0, 12.0, 100)
    mag_hb = round.(cef_magnetization(ion, blm, fields, "y", T), digits=7)
    mag_hc = round.(cef_magnetization(ion, blm, fields, "z", T), digits=7)

    f = Figure()
    Axis(f[1, 1],
        title="Yb3+ in hexagonal environment: Magnetization, " * "T=$(T)K",
        xlabel = "Field [Tesla]",
        ylabel = "Magnetic moment [Bohr magneton]"
        )
    lines!(f[1, 1], fields, mag_hb[:, 2], label=L"\mu_{0}H\perp c")
    lines!(f[1, 1], fields, mag_hc[:, 3], label=L"\mu_{0}H\parallel c")
    axislegend()
    xlims!(0, 12)
    ylims!(0, 3.5)
    f
end


"""
reproduces figure 14, page 222
"""
function plot_chi(ion, blm)
    mult_fact = 1.0
    temps = LinRange(0.5, 300, 100)
    bx = 0; by = 0.05; bz = 0
    chi_hb = cef_susceptibility(ion, blm, temps, bx, by, bz) .*mult_fact
    bx = 0; by = 0; bz = 0.05
    chi_hc = cef_susceptibility(ion, blm, temps, bx, by, bz) .*mult_fact

    f = Figure()
    Axis(f[1, 1],
        title="Yb3+ in hexagonal environment: Susceptibility, " * "B=0.05 T",
        xlabel = "Temperature [Kelvin]",
        ylabel = "Inverse mag. moment [Bohr magneton]"
        )
    lines!(f[1, 1], temps, 1 ./chi_hb[2, 2, :], label=L"\mu_{0}H\perp c")
    lines!(f[1, 1], temps, 1 ./chi_hc[3, 3, :], label=L"\mu_{0}H\parallel c")
    lines!(f[1, 1], temps, 1 ./((2*chi_hb[2, 2, :] + chi_hc[3, 3, :])./3), label=L"\chi_{\text{total}}")
    axislegend()
    xlims!(0, 300)
    # ylims!(0, 3.5)
    f
end


"""
reproduces figure 8, page 215
"""
function plot_cv(ion, blm)
    temps = LinRange(0.5, 300, 100)
    cv = cef_heatcap(ion, blm, temps, 1:4)
    cv *= 1.8R

    f = Figure()
    Axis(f[1, 1],
        title="Yb3+ in hexagonal environment: Specific heat capacity",
        xlabel = "Temperature [Kelvin]",
        ylabel = "Schottky Heat Capacity [J/molK]"
        )
    lines!(f[1, 1], temps, cv, label="Model A")
    axislegend()
    xlims!(0, 300)
    ylims!(0, 8)
    f
end


# plot_mag(yb, bdf)

# plot_chi(yb, bdf)

plot_cv(yb, bdf)