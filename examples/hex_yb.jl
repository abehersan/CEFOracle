"""
Reproduction of results from Rotter, Bauer paper, CEF model A
"""


using CEFOracle, LaTeXStrings, Plots
gr()


println("bruh")
yb = single_ion("Yb3")
bs_modelA = Dict(
        "B20"=>0.5622, "B40"=>1.6087e-5, "B60"=>6.412e-7, "B66"=>-8.324e-6)
bs_modelB = Dict(
        "B20"=>-0.06364, "B40"=>2.713e-3, "B60"=>-5.5335e-6, "B66"=>-6.659e-5)
bdf_A = blm_dframe(bs_modelA)
bdf_B = blm_dframe(bs_modelB)


"""
calculate the inelastic neutron scattering x-section at different temperatures
with model A, reproduces figure 7 of the paper
"""
es = LinRange(-12, 12, 150)
# es = LinRange(10, 150, 700)
temps = [10.0, 50.0, 200.0]
ins_plot = plot(xlabel="Energy [meV]", ylabel="I(Q, E) [arb. units]")
for t in temps
    ins_xsection = [cef_neutronxsection(yb, bdf_A, e, 2.55, t) for e in es]
    plot!(es, ins_xsection, label="T=$t K")
end


"""
calculate the Schottky specific heat capacity for both models
reproduces figure 8 of the paper
"""
temps = LinRange(0.5, 300, 150)
# cv_A = [cef_heatcapacity_speclevels(yb,bdf_A,t,1:6) for t in temps]
# cv_B = [cef_heatcapacity_speclevels(yb,bdf_B,t,1:6) for t in temps]
cv_A = [cef_heatcapacity(yb, bdf_A, t) for t in temps]
cv_B = [cef_heatcapacity(yb, bdf_B, t) for t in temps]
cv_plot = plot(
            temps, [cv_A cv_B],
            label=["A" "B"],
            xlabel="Temperature [Kelvin]",
            ylabel = "Cv [J/K/mol]"
           )


"""
calculate the magnetic moment along an applied field's direction
parallel and perpendicular to the crystallographic c-axis
reproduces figure 11 of the paper
"""
T = 0.5
fields = LinRange(0.0, 12, 150)
mag_A_para = [cef_magnetization(yb,bdf_A,T,[b,0,0],"atomic")[1] for b in fields]
mag_A_perp = [cef_magnetization(yb,bdf_A,T,[0,0,b],"atomic")[3] for b in fields]
mag_B_para = [cef_magnetization(yb,bdf_B,T,[b,0,0],"atomic")[1] for b in fields]
mag_B_perp = [cef_magnetization(yb,bdf_B,T,[0,0,b],"atomic")[3] for b in fields]
mag_A_powd = [cef_magnetization(yb,bdf_A,T,b,"atomic") for b in fields]
# mag_A_powd = 2/3 .* mag_A_perp + 1/3 .* mag_A_perp
mag_B_powd = [cef_magnetization(yb,bdf_B,T,b,"atomic") for b in fields]
mag_plot = plot(
            fields, [mag_A_para mag_A_perp mag_B_para mag_B_perp mag_A_powd mag_B_powd],
            label=[L"B\parallel c_{A}" L"B\perp c_{A}" L"B\parallel c_{B}" L"B\perp c_{B}" "Powder A" "Powder B"],
            xlabel="Field [Tesla]",
            ylabel="M " * L"[\mu_{B}/ion]"
            )


"""
calculate the inverse susceptibility for model A
two field directions are considered, the strength of the field is 0.05 Tesla
in addition, a powder average is calculated as (2chi_aa + chi_c)/3
reproduces figure 14 of the paper
"""
bext = 0.05
temps = LinRange(0.5, 300, 150)
invchi_para = [1/cef_susceptibility(yb,bdf_A,t,[0,0,bext],"CGS")[3] for t in temps]
invchi_perp = [1/cef_susceptibility(yb,bdf_A,t,[bext,0,0],"CGS")[1] for t in temps]
invchi_powd = [1/cef_susceptibility(yb,bdf_A,t,bext, "CGS") for t in temps]
invchi_powd_c = (2/3 .* invchi_perp + 1/3 .* invchi_para)
invchi_plot = plot(
                temps, [invchi_para invchi_perp invchi_powd invchi_powd_c],
                label=[L"B\parallel c" L"B\perp c" "Powder" "Calc powder"],
                xlabel="Temperature [Kelvin]",
                ylabel=L"1/\chi\;[mol/emu]"
                )


plot(ins_plot, cv_plot, mag_plot, invchi_plot, layout=(2, 2), legend=true,
    dpi=300, size=(1080, 720), thickness_scaling=1.5, scalefontsizes=1.5)
