"""
cef_fit.jl

Extract crystal-field parameters from data (magnetization, susceptibility,
specific heat capacity and INS)
"""


using CEFOracle, DataFrames, CSV, StatsPlots, LaTeXStrings
using Optimization, OptimizationNLopt, OptimizationOptimJL


function dset_model_plot(; ion::mag_ion, blm::DataFrame, dsets::cef_datasets)
    chi_plot = plot(xlabel="Temperature [K]", ylabel="Chi [cm^3/mol]")
    mag_plot = plot(xlabel="Field [T]", ylabel="Mag. mom. [BM/f.u.]")
    # TODO: iterate over the groups in each entry in the datasets
    # struct, --- find a way of iterating over the structs' items more
    # easily!

    # if !isnothing(dsets.mag_data)
    #     @df dsets.mag_data scatter!(mag_plot, :Bext, :MB, yerr=:Err, group=:Dir,
    #                                 m=:x, la=0)
    #     temps = Set(dsets.mag_data.T)
    #     for t in temps
    #         plot!()
    #     end
    # end
    # if !isnothing(dsets.susc_data)
    #     @df dsets.susc_data scatter!(chi_plot, :T, :MChi, yerr=:Err,
    #                                  group=:Dir, m=:x, la=0)
    # end
    l = @layout [a ; c]
    plots = [chi_plot, mag_plot]
    xlims!(chi_plot, 0, 25)
    fit_plot = plot(plots..., layout=l, legend=:bottomright, dpi=300,
                    size=(1080,720), thickness_scaling=1.5, scalefontsizes=1.5)
    fit_plot
end


function callback(p, chi2)
    false
end


function cef_objective(u::Vector, p::Vector)::Real
    # u = vector of stevens parameters in the same order as blm_symbols
    # p = [ion, cef_data, blm_symbols, np]
    ion::mag_ion, cef_data::cef_datasets, blm_symbols::Vector, np::Int = p
    bs_dframe = blm_dframe(Dict((s, v) for (s, v) in zip(blm_symbols, u)))
    chi2_cef(ion, bs_dframe, cef_data)/(np - length(u))
    # chi2_cef(ion, bs_dframe, cef_data)
end


ho = single_ion("Ho3")
initial_bs = Dict("B20"=>-0.05575208878804679,
                  "B40"=>-4.636243790111825e-6,
                  "B43"=> -8.660454382441345e-7,
                  "B4m3"=>-3.870452739528808e-7,
                  "B60"=>-0.0025412414866357196,
                  "B63"=>3.087773150542659e-5,
                  "B6m3"=>-0.00011623123127658234,
                  "B66"=>-0.011212321320726323,
                  "B6m6"=>8.75185498259519e-8)
bdf = blm_dframe(initial_bs)
b_pars = collect(values(initial_bs))
b_syms = collect(keys(initial_bs))
b_ubounds = [3.5 for b in b_pars]
b_lbounds = [-3.5 for b in b_pars]


data_path::String = "../../../DataAnalysis/CEF/artifacts/ext_data/"
hobr3_chi_data = DataFrame(CSV.File(data_path*"HoBr3/HoBr3_chi_all.csv"))
hobr3_mag_data = DataFrame(CSV.File(data_path*"HoBr3/HoBr3_mag_all.csv"))
hobr3_mpms_data = cef_datasets(mag_data=hobr3_mag_data,
                               susc_data=hobr3_chi_data)
np = nrow(hobr3_chi_data) + nrow(hobr3_mag_data)
p = [ho, hobr3_mpms_data, b_syms, np]


println("Initial Chi2/NDOF: $(cef_objective(b_pars, p))")
dset_model_plot(dsets=hobr3_mpms_data)


# opt_func = OptimizationFunction(cef_objective, syms=b_syms)
# opt_prob = OptimizationProblem(opt_func, b_pars, p, sense=:MinSense,
                               # lb=b_lbounds, ub=b_ubounds)
# opt_sol = solve(opt_prob, Optim.NelderMead(), show_trace=true, show_every=25)
# opt_sol = solve(opt_prob, Optim.ParticleSwarm(lower=b_lbounds, upper=b_ubounds,
#                                              n_particles=70),
#                show_trace=true, show_every=15, maxiters=1e6, maxtime=Inf)
