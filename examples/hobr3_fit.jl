"""
cef_fit.jl

Extract crystal-field parameters from data (magnetization and susceptibility)
Example of HoBr3
"""


using CEFOracle, DataFrames, CSV, StatsPlots, LaTeXStrings
using Optimization, OptimizationNLopt, OptimizationOptimJL


function plot_chi(cef_data::cef_datasets)
    chi_data::DataFrame = cef_data.chi_data
    gchi_data = groupby(chi_data, :Dir)
    ion::mag_ion = cef_data.ion; blm::DataFrame = cef_data.blm
    chi_plot = plot(xlabel="Temperature [K]", ylabel="Chi [cm^3/mol]")
    @df chi_data scatter!(chi_plot, :Chi, :T, yerr=:Err, group=:Dir, m=:x)
    for ((; dir), single_df) in pairs(gchi_data)
        x_calc = single_df.T
        y_calc = zeros(nrow(single_df))
        for (i, pnt) in enumerate(eachrow(single_df))
            bext = begin
                if dir == "a"
                    [pnt.B, 0, 0]
                elseif dir == "b"
                    [0, pnt.B, 0]
                elseif dir == "c"
                    [0, 0, pnt.B]
                elseif dir == "p"
                    pnt.B
                end
            end
            y_calc[i] = cef_susceptibility(ion, blm, pnt.T, bext)
        end
        plot!(chi_plot, x_calc, y_calc, label="$dir")
    end
    chi_plot
end


function callback(p, chi2, cef_data)
    false
end


function cef_objective(u::Vector, p::Vector)::Real
    # u = vector of stevens parameters in the same order as blm_symbols
    # p = [ion, cef_data, blm_symbols, np]
    ion::mag_ion, cef_data::cef_datasets, blm_symbols::Vector, np::Int = p
    bs_dframe = blm_dframe(Dict((s, v) for (s, v) in zip(blm_symbols, u)))
    chi2_cef(cef_data)/(np - length(u)), cef_data
end


# function main()
    ho = single_ion("Ho3")
    initial_bs = Dict("B20"=> -0.05575208878804679,
                      "B40"=> -4.636243790111825e-6,
                      "B43"=> -8.660454382441345e-7,
                      "B4m3"=>-3.870452739528808e-7,
                      "B60"=> -0.0025412414866357196,
                      "B63"=> +3.087773150542659e-5,
                      "B6m3"=>-0.00011623123127658234,
                      "B66"=> -0.011212321320726323,
                      "B6m6"=>+8.75185498259519e-8)
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
    p = [ho, hobr3_all_data, b_syms, np]

    println("Initial Chi2/NDOF: $(cef_objective(b_pars, p))[1]")
    plot_chi(hobr3_all_data)

    # opt_func = OptimizationFunction(cef_objective, syms=b_syms)
    # opt_prob = OptimizationProblem(opt_func, b_pars, p, sense=:MinSense,
                                # lb=b_lbounds, ub=b_ubounds)
    # opt_sol = solve(opt_prob, Optim.NelderMead(), show_trace=true, show_every=25)
    # opt_sol = solve(opt_prob, Optim.ParticleSwarm(lower=b_lbounds, upper=b_ubounds,
    #                                              n_particles=70),
    #                show_trace=true, show_every=15, maxiters=1e6, maxtime=Inf)

# end


# main()