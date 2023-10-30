function free_energy(Ep::Vector{Float64}, T::Real)::Float64
    return kB*T*partition_function(Ep, T)
end


function calc_heatcap(; Ep::Vector{Float64}, T::Real)::Float64
    np = population_factor(Ep, T)
    heatcap = sum((Ep/(kB*T)).^2 .* np)
    heatcap -= sum((Ep/(kB*T).*np).^2)
    return heatcap
end


function calc_entropy(; Ep::Vector{Float64}, T::Real)::Float64
    Z = partition_function(Ep, T)
    return (log2(Z) - log2(partition_function(Ep, 1e-3)))
end


function cef_entropy!(single_ion::mag_ion, bfactors::DataFrame,
                     calc_grid::DataFrame; units::String="SI")::Nothing
    convfac = begin
        if isequal(units, "SI")
            NA * 1.602176487*1e-22 *kB  # NA * muB [J/mol/K]
        else
            @error "Units not supported in Cv calculations. Use SI, [J/K/mol]"
        end
    end
    E = eigvals(cef_hamiltonian(single_ion, bfactors))
    E .-= minimum(E)
    calc_grid[!, :CP_CALC] = [calc_heatcap(Ep=E, T=t)*convfac for t in calc_grid.T]
    calc_grid[!, :S_CALC]  = [calc_entropy(Ep=E, T=t)*convfac for t in calc_grid.T]
    return
end


function cef_entropy_speclevels!(single_ion::mag_ion, bfactors::DataFrame,
                                calc_grid::DataFrame; levels::UnitRange=1:4,
                                units::String="SI")::Nothing
    # only levels specified contribute (2J+1 levels total)
    convfac = begin
        if isequal(units, "SI")
            NA * 1.602176487*1e-22 *kB  # NA * muB [J/mol/K]
        else
            @error "Units not supported in Cv calculations. Use SI, [J/K/mol]"
        end
    end
    E = eigvals(cef_hamiltonian(single_ion, bfactors))
    E .-= minimum(E)
    calc_grid[!, :CP_CALC] = [calc_heatcap(Ep=E[levels], T=t)*convfac for t in calc_grid.T]
    calc_grid[!, :S_CALC]  = [calc_entropy(Ep=E[levels], T=t)*convfac for t in calc_grid.T]
    return
end