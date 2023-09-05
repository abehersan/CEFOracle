"""
Reproduction of results from the paper by Rotter and Bauer
"""


using CEFOracle, LaTeXStrings, Plots


yb = single_ion("Yb3")
bfactors_A = blm_dframe(Dict("B20"=>0.5622, "B40"=>1.6087e-5, "B60"=>6.412e-7, "B66"=>-8.324e-6))
bfactors_B = blm_dframe(Dict("B20"=>-0.06364, "B40"=>2.713e-3, "B60"=>-5.5335e-6, "B66"=>-6.659e-5))


"""
calculate the inelastic neutron scattering x-section at different temperatures
with model A, reproduces figure 7 of the paper
"""
es = LinRange(-12, 12, 700)
temps = [10.0, 50.0, 200.0]
ins_plot = plot(xlabel="Energy [meV]", ylabel="I(Q, E) [arb. units]")
for t in temps
    ins_xsection = zeros(Float64, length(es))
    for i in eachindex(es)
        ins_xsection[i] = cef_neutronxsection_powder(yb, bfactors_A, E=es[i],
                                                    Q=2.25, T=t)
    end
    plot!(es, ins_xsection, label="T=$t K")
end


"""
calculate the Schottky specific heat capacity for both models
reproduces figure 8 of the paper
"""
temps = LinRange(0.5, 300, 150)
cv_A = zeros(Float64, length(temps))
cv_B = zeros(Float64, length(temps))
for i in eachindex(temps)
    cv_A[i] = cef_heatcapacity(yb, bfactors_A, T=temps[i])
    cv_B[i] = cef_heatcapacity(yb, bfactors_B, T=temps[i])
end
cv_plot = plot(temps, [cv_A cv_B], label=["A" "B"],
               xlabel="Temperature [Kelvin]", ylabel = "Cv [J/K/mol]")


"""
calculate the magnetic moment along an applied field's direction
parallel and perpendicular to the crystallographic c-axis
reproduces figure 11 of the paper
"""
T = 0.5
fields = LinRange(0.0, 12, 150)
mag_A_para = zeros(Float64, length(fields))
mag_A_perp = zeros(Float64, length(fields))
mag_B_para = zeros(Float64, length(fields))
mag_B_perp = zeros(Float64, length(fields))
for i in eachindex(fields)
    mag_A_para[i] = cef_magnetization_crystal(yb, bfactors_A, T=T, B=[fields[i],0,0], units="ATOMIC")
    mag_A_perp[i] = cef_magnetization_crystal(yb, bfactors_A, T=T, B=[0,0,fields[i]], units="ATOMIC")
    mag_B_para[i] = cef_magnetization_crystal(yb, bfactors_B, T=T, B=[fields[i],0,0], units="ATOMIC")
    mag_B_perp[i] = cef_magnetization_crystal(yb, bfactors_B, T=T, B=[0,0,fields[i]], units="ATOMIC")
end
mag_plot = plot(fields, [mag_A_para mag_A_perp mag_B_para mag_B_perp],
            label=[L"B\parallel c_{A}" L"B\perp c_{A}" L"B\parallel c_{B}" L"B\perp c_{B}"],
            xlabel="Field [Tesla]", ylabel="M " * L"[\mu_{B}/ion]")


"""
calculate the inverse susceptibility for model A
two field directions are considered, the strength of the field is 0.05 Tesla
in addition, a powder average is calculated as (2chi_aa + chi_c)/3
reproduces figure 14 of the paper
"""
bext = 0.05
temps = LinRange(0.5, 300, 150)
invchi_para = zeros(Float64, length(temps))
invchi_perp = zeros(Float64, length(temps))
for i in eachindex(temps)
    invchi_para[i] = 1/cef_susceptibility_crystal(yb, bfactors_A, T=temps[i], B=[0,0,bext], units="CGS")
    invchi_perp[i] = 1/cef_susceptibility_crystal(yb, bfactors_A, T=temps[i], B=[bext,0,0], units="CGS")
end
invchi_powd_c = (2/3 .* invchi_perp .+ 1/3 .* invchi_para)
invchi_plot = plot(temps, [invchi_para invchi_perp invchi_powd_c],
                label=["B parallel c" "B perp c" "Powder"],
                xlabel="Temperature [Kelvin]", ylabel="1/Susceptibility [mol/emu]")


plot(ins_plot, cv_plot, mag_plot, invchi_plot, layout=(2, 2), legend=true,
    dpi=300, size=(1080, 720), thickness_scaling=1.5, scalefontsizes=1.5)
