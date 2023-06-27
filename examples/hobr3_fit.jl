"""
cef_fit.jl

Extract crystal-field parameters from data (magnetization and susceptibility)
Example of HoBr3
"""


using CEFOracle, DataFrames, CSV, StatsPlots, LaTeXStrings
using Optimization, OptimizationNLopt, OptimizationOptimJL


function plot_chi(cef_data::cef_datasets)
    chi_data::DataFrame = cef_data.chi_data
    ion::mag_ion = cef_data.ion; blm::DataFrame = cef_data.blm
    chi_plot = plot(xlabel="Temperature [K]", ylabel="Chi [cm^3/mol]")
    @df chi_data scatter!(chi_plot, :T, :Chi, yerr=:Err, group=:Dir, m=:x)
    for ((; Dir), single_df) in pairs(groupby(chi_data, :Dir))
        x_calc = single_df.T
        y_calc = zeros(nrow(single_df))
        for (i, pnt) in enumerate(eachrow(single_df))
            y_calc[i] = begin
                if Dir == "a"
                    cef_susceptibility(ion,blm,pnt.T,[pnt.Bext,0,0],"CGS")[1]
                elseif Dir == "b"
                    cef_susceptibility(ion,blm,pnt.T,[0,pnt.Bext,0],"CGS")[2]
                elseif Dir == "c"
                    cef_susceptibility(ion,blm,pnt.T,[0,0,pnt.Bext],"CGS")[3]
                elseif Dir == "p"
                    cef_susceptibility(ion,blm,pnt.T,pnt.Bext,"CGS")
                end
            end
        end
        plot!(chi_plot, x_calc, y_calc, label="$Dir")
    end
    display(chi_plot)
end


function callback(p, chi2, cef_data)
    false
end


function cef_objective(u::Vector, p::Tuple)::Real
    # u = vector of stevens parameters in the same order as blm_symbols
    # p = [ion, cef_data, blm_symbols, np]
    ion::mag_ion, cef_data::cef_datasets, blm_symbols::Vector, np::Int = p
    bs_dframe = blm_dframe(Dict((s, v) for (s, v) in zip(blm_symbols, u)))
    chi2_cef(cef_data)/(np - length(u))#, cef_data
    # chi2_cef(cef_data)
end


ho = single_ion("Ho3")
initial_bs = Dict("B20"=>-0.007,
                 "B4m3"=>-1.00256986,
                 "B40"=>-6.45614e-05,
                 "B43"=>5.45495e-06,
                 "B6m6"=>5.9241e-09,
                 "B6m3"=>-1.30549e-02,
                 "B60"=>1.47946e-07,
                 "B63"=>1.24976e-08,
                 "B66"=>-1.0817e-06,)
bdf = blm_dframe(initial_bs)
b_pars = collect(values(initial_bs))
b_syms = collect(keys(initial_bs))
b_ubounds = [3.5 for b in b_pars]
b_lbounds = [-3.5 for b in b_pars]

data_path::String = "../../../DataAnalysis/CEF/artifacts/ext_data/"
hobr3_chi_data = DataFrame(CSV.File(data_path*"HoBr3/HoBr3_chi_all.csv"))
hobr3_mag_data = DataFrame(CSV.File(data_path*"HoBr3/HoBr3_mag_all.csv"))
hobr3_all_data = cef_datasets(ion=ho, blm=bdf, mag_data=hobr3_mag_data,
                              chi_data=hobr3_chi_data)
np = nrow(hobr3_chi_data) + nrow(hobr3_mag_data)
p = (ho, hobr3_all_data, b_syms, np)

println("Initial Chi2/NDOF: $(cef_objective(b_pars, p))")
# plot_chi(hobr3_all_data)

opt_func = OptimizationFunction(cef_objective, syms=b_syms)
opt_prob = OptimizationProblem(opt_func, b_pars, p, sense=:MinSense,)
#                             # lb=b_lbounds, ub=b_ubounds)
# opt_sol = solve(opt_prob, Optim.NelderMead(), show_trace=true, show_every=25)
opt_sol = solve(opt_prob, Optim.ParticleSwarm(lower=b_lbounds, upper=b_ubounds,
                                             n_particles=50),
               show_trace=true, show_every=15, maxiters=1e6, maxtime=Inf)